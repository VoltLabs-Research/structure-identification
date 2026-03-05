#include <volt/core/volt.h>
#include <volt/utilities/concurrence/parallel_system.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/polyhedral_template_matching.h>
#include <volt/analysis/ptm_neighbor_finder.h>
#include <ptm_constants.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <numeric>

namespace Volt {

StructureAnalysis::StructureAnalysis(
    AnalysisContext& context,
    bool identifyPlanarDefects, 
    Mode identificationMode,
    float rmsd
) :
    _context(context),
    _identificationMode(identificationMode),
    _rmsd(rmsd),    
    _clusterGraph(std::make_unique<ClusterGraph>()),
    _coordStructures(
        _context.structureTypes, 
        context.inputCrystalType, 
        identifyPlanarDefects, 
        context.simCell
    )
{
    static std::once_flag init_flag;
    std::call_once(init_flag, []() {
        CoordinationStructures::initializeStructures();
    });

    int requestedMaxNeighbors = 0;
    if(usingPTM()){
        requestedMaxNeighbors = PTM_MAX_NBRS;
    }else{
        requestedMaxNeighbors = std::max(
            _coordStructures.getLatticeStruct(_context.inputCrystalType).maxNeighbors,
            _coordStructures.getCoordinationNumber()
        );
        if(requestedMaxNeighbors <= 0) requestedMaxNeighbors = 1;
    }

    _context.atomSymmetryPermutations = std::make_shared<ParticleProperty>(
        _context.atomCount(), DataType::Int, 1, 0, true);

    _context.neighborOffsets = std::make_shared<ParticleProperty>(
        _context.atomCount() + 1, DataType::Int, 1, 0, true);

    // O(1) neighbor count per atom — avoids linear -1 sentinel scan
    _context.neighborCounts = std::make_shared<ParticleProperty>(
        _context.atomCount(), DataType::Int, 1, 0, true);

    _context.templateIndex = std::make_shared<ParticleProperty>(
        _context.atomCount(), DataType::Int, 1, 0, true);

    std::fill(_context.neighborOffsets->dataInt(),
              _context.neighborOffsets->dataInt() + _context.neighborOffsets->size(),
              0);
    std::fill(_context.structureTypes->dataInt(),
              _context.structureTypes->dataInt() + _context.structureTypes->size(),
              LATTICE_OTHER);

}

json StructureAnalysis::getPerAtomProperties(
    const LammpsParser::Frame &frame,
    const std::vector<int>* structureTypes
){
    json perAtom = json::array();

    for(size_t i = 0; i < frame.natoms; ++i){
        int structureType = 0;
        if(structureTypes && i < structureTypes->size()){
            structureType = (*structureTypes)[i];
        }

        json atom;
        atom["id"] = i < frame.ids.size() ? frame.ids[i] : static_cast<int>(i);
        atom["structure_type"] = structureType;
        atom["structure_name"] = getStructureTypeName(structureType);

        if(i < static_cast<size_t>(frame.positions.size())){
            const auto &pos = frame.positions[i];
            atom["pos"] = {pos.x(), pos.y(), pos.z()};
        }else{
            atom["pos"] = {0.0, 0.0, 0.0};
        }

        perAtom.push_back(atom);
    }

    return perAtom;
}

bool StructureAnalysis::setupPTM(Volt::PTM& ptm, size_t N){
	ptm.setCalculateDefGradient(false);
	ptm.setRmsdCutoff(std::numeric_limits<double>::infinity());
	
	return ptm.prepare(_context.positions->constDataPoint3(), N, _context.simCell);
}

void StructureAnalysis::storeNeighborIndices(const PTM::Kernel& kernel, size_t atomIndex){
    int numNeighbors = kernel.numTemplateNeighbors();
    _context.neighborCounts->setInt(atomIndex, numNeighbors);
}

void StructureAnalysis::storeOrientationData(const PTM::Kernel& kernel, size_t atomIndex){
    auto quaternion = kernel.orientation();
    double* orientation = _context.ptmOrientation->dataDouble() + 4 * atomIndex;
    
    orientation[0] = quaternion.x();
    orientation[1] = quaternion.y();
    orientation[2] = quaternion.z();
    orientation[3] = quaternion.w();
}

void StructureAnalysis::storeDeformationGradient(const PTM::Kernel& kernel, size_t atomIndex) {
	if(!_context.ptmDeformationGradient){
		return;
	}
	const auto& F = kernel.deformationGradient();
	double* F_dest = _context.ptmDeformationGradient->dataDouble() + 9 * atomIndex;
    const double* F_src = F.elements();
    
    for(int k = 0; k < 9; ++k){
        F_dest[k] = F_src[k];
    }
}

int StructureAnalysis::findClosestSymmetryPermutation(int structureType, const Matrix3& rotation){
    const LatticeStructure& lattice = CoordinationStructures::getLatticeStruct(structureType);
    int bestIndex = 0;
    double bestDeviation = std::numeric_limits<double>::max();

    for(int i = 0; i < lattice.permutations.size(); ++i){
        const Matrix3& sym = lattice.permutations[i].transformation;
        double deviation = 0;
        for(int r = 0; r < 3; ++r){
            for(int c = 0; c < 3; ++c){
                double diff = rotation(r, c) - sym(r, c);
                deviation += diff * diff;
            }
        }
        if(deviation < bestDeviation){
            bestDeviation = deviation;
            bestIndex = i;
        }
    }
    return bestIndex;
}

// Compute the maximum distance of any neighbor from a crystalline atom
void StructureAnalysis::computeMaximumNeighborDistanceFromPTM(){
    const size_t N = _context.atomCount();
    if(N == 0){
        _maximumNeighborDistance = 0.0;
        return;
    }

    const auto* pos = _context.positions->constDataPoint3();
    const auto& invMat = _context.simCell.inverseMatrix();
    const auto& dirMat = _context.simCell.matrix();
    const int* counts = _context.neighborCounts->constDataInt();
    const int* offsets = _context.neighborOffsets->constDataInt();
    const int* indices = _context.neighborIndices->constDataInt();

    double maxDistance = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, N),
        0.0,
        [&](const tbb::blocked_range<size_t>& r, double max_so_far) -> double {
            for(size_t i = r.begin(); i < r.end(); ++i){
                const int nNeigh = counts[i];
                double localMaxDist = 0.0;
                const int start = offsets[i];
                for(int j = 0; j < nNeigh; ++j){
                    int nb = indices[start + j];

                    Vector3 delta = pos[nb] - pos[i];
                    double f[3];
                    for(int d = 0; d < 3; ++d){
                        f[d] = invMat.prodrow(delta, d);
                        f[d] -= std::round(f[d]);
                    }

                    Vector3 mind;
                    mind = dirMat.column(0) * f[0];
                    mind += dirMat.column(1) * f[1];
                    mind += dirMat.column(2) * f[2];
                    double d = mind.length();
                    if(d > localMaxDist) localMaxDist = d;
                }
                if(localMaxDist > max_so_far) max_so_far = localMaxDist;
            }
            return max_so_far;
        },
        [](double a, double b) -> double { return std::max(a, b); }
    );

    spdlog::debug("Maximum neighbor distance (from PTM): {}", maxDistance);
    _maximumNeighborDistance = maxDistance;
}

// Runs the Polyhedral Template Matching (PTM) algorithm on every atom,
// collects raw RMSD values (with no initial cutoff).
void StructureAnalysis::determineLocalStructuresWithPTM() {
    const size_t N = _context.atomCount();
    if(!N) return;

    Volt::PTM ptm;
    if(!setupPTM(ptm, N)){
        throw std::runtime_error("Error trying to initialize PTM.");
    }

    _context.ptmOrientation = std::make_shared<ParticleProperty>(N, DataType::Double, 4, 0.0, true);
    _context.ptmRmsd = std::make_shared<ParticleProperty>(N, DataType::Double, 1, 0.0, true);
    _context.correspondencesCode = std::make_shared<ParticleProperty>(N, DataType::Int64, 1, 0, true);
    if(ptm.calculateDefGradient()){
        _context.ptmDeformationGradient = std::make_shared<ParticleProperty>(N, DataType::Double, 9, 0.0, true);
    }else{
        _context.ptmDeformationGradient.reset();
    }

    // Clear arrays for second pass
    std::fill(_context.neighborCounts->dataInt(),
              _context.neighborCounts->dataInt() + _context.neighborCounts->size(), 0);
    std::fill(_context.structureTypes->dataInt(),
              _context.structureTypes->dataInt() + _context.structureTypes->size(), LATTICE_OTHER);


    std::vector<uint64_t> cached(N, 0ull);
    std::vector<int> localCounts(N, 0);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto &r){
        PTM::Kernel kernel(ptm);
        for(size_t i = r.begin(); i < r.end(); ++i){
            kernel.cacheNeighbors(i, &cached[i]);
            StructureType type = kernel.identifyStructure(i, cached);

            const double rmsd = kernel.rmsd();
            _context.ptmRmsd->setDouble(i, rmsd);

            auto* c = reinterpret_cast<uint64_t*>(_context.correspondencesCode->data());
            c[i] = kernel.correspondencesCode();

            // Only keep atoms whose RMSD <= finalCutoff
            if(type == StructureType::OTHER || rmsd > _rmsd){
                continue;
            }

            _context.structureTypes->setInt(i, type);
            const int numNeighbors = kernel.numTemplateNeighbors();
            localCounts[i] = numNeighbors;
            _context.neighborCounts->setInt(i, numNeighbors);
            storeOrientationData(kernel, i);
            storeDeformationGradient(kernel, i);
            _context.templateIndex->setInt(i, kernel.bestTemplateIndex());
        }
    });

    auto* offsets = _context.neighborOffsets->dataInt();
    offsets[0] = 0;
    for(size_t i = 0; i < N; ++i){
        offsets[i + 1] = offsets[i] + localCounts[i];
    }

    const size_t totalNeighbors = static_cast<size_t>(offsets[N]);
    _context.neighborIndices = std::make_shared<ParticleProperty>(
        totalNeighbors, DataType::Int, 1, 0, false);
    auto* indices = _context.neighborIndices->dataInt();

    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto& r){
        PTM::Kernel kernel(ptm);
        for(size_t i = r.begin(); i < r.end(); ++i){
            const int count = localCounts[i];
            if(count == 0){
                continue;
            }
            kernel.cacheNeighbors(i, &cached[i]);
            kernel.identifyStructure(i, cached);
            const int start = offsets[i];
            for(int j = 0; j < count; ++j){
                indices[start + j] = kernel.getTemplateNeighbor(j).index;
            }
        }
    });

    for(size_t i = 0; i < N; ++i){
        if(_context.neighborCounts->getInt(i) == 0){
            _context.neighborCounts->setInt(i, localCounts[i]);
        }
    }
}

void StructureAnalysis::identifyStructuresCNA(){
    const int maxNeighborListSize = MAX_NEIGHBORS;
    NearestNeighborFinder neighFinder(maxNeighborListSize);
    if(!neighFinder.prepare(_context.positions, _context.simCell, _context.particleSelection)){
        throw std::runtime_error("Error in neighFinder.preapre(...)");
    }

    const size_t N = _context.atomCount();
    std::vector<int> localCounts(N, 0);
    _maximumNeighborDistance = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, N),
        0.0, [this, &neighFinder, &localCounts](const tbb::blocked_range<size_t>& r, double max_dist_so_far) -> double {
            for(size_t index = r.begin(); index != r.end(); ++index){
                int count = 0;
                double localMaxDistance = _coordStructures.determineLocalStructure(neighFinder, index, &count);
                localCounts[index] = count;
                if (localMaxDistance > max_dist_so_far) {
                    max_dist_so_far = localMaxDistance;
                }
            }
            return max_dist_so_far;
        },
        [](double a, double b) -> double {
            return std::max(a, b);
        }
    );

    auto* offsets = _context.neighborOffsets->dataInt();
    offsets[0] = 0;
    for(size_t i = 0; i < N; ++i){
        offsets[i + 1] = offsets[i] + localCounts[i];
    }
    const size_t totalNeighbors = static_cast<size_t>(offsets[N]);
    _context.neighborIndices = std::make_shared<ParticleProperty>(
        totalNeighbors, DataType::Int, 1, 0, false);

    auto* indices = _context.neighborIndices->dataInt();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto& r){
        NearestNeighborFinder::Query<MAX_NEIGHBORS> query(neighFinder);
        for(size_t index = r.begin(); index != r.end(); ++index){
            const int count = localCounts[index];
            if(count <= 0){
                continue;
            }
            query.findNeighbors(neighFinder.particlePos(index));
            const int start = offsets[index];
            for(int j = 0; j < count; ++j){
                indices[start + j] = query.results()[j].index;
            }
            _context.neighborCounts->setInt(index, count);
        }
    });
}

void StructureAnalysis::appendNeighbors(const std::vector<std::vector<int>>& extras){
    if(extras.empty()){
        return;
    }

    const size_t N = _context.atomCount();
    if(extras.size() != N){
        throw std::runtime_error("appendNeighbors: size mismatch");
    }

    const int* oldOffsets = _context.neighborOffsets->constDataInt();
    const int* oldIndices = _context.neighborIndices->constDataInt();
    std::vector<int> newOffsets(N + 1, 0);
    for(size_t i = 0; i < N; ++i){
        const int oldCount = _context.neighborCounts->getInt(static_cast<int>(i));
        const int extraCount = static_cast<int>(extras[i].size());
        if(oldCount + extraCount > MAX_NEIGHBORS){
            throw std::runtime_error("appendNeighbors: neighbor overflow");
        }
        newOffsets[i + 1] = newOffsets[i] + oldCount + extraCount;
    }

    const size_t totalNeighbors = static_cast<size_t>(newOffsets[N]);
    auto newIndicesProp = std::make_shared<ParticleProperty>(
        totalNeighbors, DataType::Int, 1, 0, false);
    int* newIndices = newIndicesProp->dataInt();

    for(size_t i = 0; i < N; ++i){
        const int oldCount = _context.neighborCounts->getInt(static_cast<int>(i));
        const int startOld = oldOffsets[i];
        const int startNew = newOffsets[i];
        std::copy(oldIndices + startOld, oldIndices + startOld + oldCount, newIndices + startNew);
        if(!extras[i].empty()){
            std::copy(extras[i].begin(), extras[i].end(), newIndices + startNew + oldCount);
        }
        _context.neighborCounts->setInt(static_cast<int>(i), oldCount + static_cast<int>(extras[i].size()));
    }

    _context.neighborOffsets = std::make_shared<ParticleProperty>(
        N + 1, DataType::Int, 1, 0, false);
    std::copy(newOffsets.begin(), newOffsets.end(), _context.neighborOffsets->dataInt());
    _context.neighborIndices = std::move(newIndicesProp);
}

void StructureAnalysis::identifyStructures(){
    if(usingPTM()){
        determineLocalStructuresWithPTM();
        computeMaximumNeighborDistanceFromPTM();
        return;
    }

    identifyStructuresCNA();
}

}
