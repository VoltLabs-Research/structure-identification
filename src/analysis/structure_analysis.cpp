#include <volt/core/volt.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/cluster_hierarchy_rebuilder.h>
#include <spdlog/spdlog.h>

#include <algorithm>

namespace Volt {

namespace StructureAnalysisDetail{

const StructureAnalysisCrystalInfo& requireCrystalInfo(
    const std::shared_ptr<const StructureAnalysisCrystalInfo>& crystalInfoProvider
){
    if(!crystalInfoProvider){
        throw std::runtime_error("StructureAnalysis crystal information provider has not been configured.");
    }
    return *crystalInfoProvider;
}

AnalysisContext& requireAnalysisContext(StructureContext& context){
    auto* analysisContext = dynamic_cast<AnalysisContext*>(&context);
    if(!analysisContext){
        throw std::runtime_error("StructureAnalysis operation requires AnalysisContext.");
    }
    return *analysisContext;
}

const AnalysisContext& requireAnalysisContext(const StructureContext& context){
    auto* analysisContext = dynamic_cast<const AnalysisContext*>(&context);
    if(!analysisContext){
        throw std::runtime_error("StructureAnalysis operation requires AnalysisContext.");
    }
    return *analysisContext;
}

}

using namespace StructureAnalysisDetail;

StructureAnalysis::StructureAnalysis(
    StructureContext& context
) :
    _context(context),
    _clusterGraph(std::make_unique<ClusterGraph>())
{
    if(!_context.atomClusters){
        _context.atomClusters = std::make_shared<ParticleProperty>(_context.atomCount(), DataType::Int, 1, 0, true);
    }
    if(!_context.neighborOffsets){
        _context.neighborOffsets = std::make_shared<ParticleProperty>(_context.atomCount() + 1, DataType::Int, 1, 0, true);
        std::fill(_context.neighborOffsets->dataInt(), _context.neighborOffsets->dataInt() + _context.neighborOffsets->size(), 0);
    }
    if(!_context.neighborCounts){
        _context.neighborCounts = std::make_shared<ParticleProperty>(_context.atomCount(), DataType::Int, 1, 0, true);
    }
    if(!_context.atomAllowedSymmetryMasks){
        _context.atomAllowedSymmetryMasks = std::make_shared<ParticleProperty>(_context.atomCount(), DataType::Int64, 1, 0, true);
    }

    if(_context.structureTypes){
        if(!_context.atomSymmetryPermutations){
            _context.atomSymmetryPermutations = std::make_shared<ParticleProperty>(_context.atomCount(), DataType::Int, 1, 0, true);
            std::fill(
                _context.atomSymmetryPermutations->dataInt(),
                _context.atomSymmetryPermutations->dataInt() + _context.atomSymmetryPermutations->size(),
                -1
            );
        }
        std::fill(_context.structureTypes->dataInt(), _context.structureTypes->dataInt() + _context.structureTypes->size(), LATTICE_OTHER);
    }
}

Cluster* StructureAnalysis::atomCluster(int atomIndex) const{
    Cluster* cluster = clusterGraph().findCluster(_context.atomClusters->getInt(atomIndex));
    if(!cluster){
        return nullptr;
    }
    return ClusterHierarchyUtils::getParentGrain(const_cast<StructureAnalysis&>(*this), cluster);
}

int StructureAnalysis::getNeighbor(int centralAtomIndex, int neighborListIndex) const{
    assert(_context.neighborOffsets && _context.neighborIndices);
    const int count = _context.neighborCounts->getInt(centralAtomIndex);

    if(neighborListIndex < 0 || neighborListIndex >= count){
        return -1;
    }
    
    const int* offsets = _context.neighborOffsets->constDataInt();
    const int* indices = _context.neighborIndices->constDataInt();
    const int start = offsets[centralAtomIndex];

    return indices[start + neighborListIndex];
}

int StructureAnalysis::findNeighbor(int centralAtomIndex, int neighborAtomIndex) const{
    assert(_context.neighborOffsets && _context.neighborIndices);
    const int count = _context.neighborCounts->getInt(centralAtomIndex);
    const int* offsets = _context.neighborOffsets->constDataInt();
    const int* indices = _context.neighborIndices->constDataInt();
    const int start = offsets[centralAtomIndex];

    for(int index = 0; index < count; index++){
        if(indices[start + index] == neighborAtomIndex){
            return index;
        }
    }

    return -1;
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
    std::vector<std::vector<int>> filteredExtras(N);
    std::size_t trimmedNeighborCount = 0;
    std::size_t affectedAtomCount = 0;

    std::vector<int> newOffsets(N + 1, 0);
    for(size_t i = 0; i < N; ++i){
        const int oldCount = _context.neighborCounts->getInt(static_cast<int>(i));
        if(oldCount > MAX_NEIGHBORS){
            throw std::runtime_error("appendNeighbors: base neighbor overflow");
        }

        const int startOld = oldOffsets[i];
        auto& filtered = filteredExtras[i];
        filtered.reserve(extras[i].size());

        for(int neighbor : extras[i]){
            if(neighbor < 0){
                continue;
            }

            bool alreadyPresent = false;
            for(int slot = 0; slot < oldCount; ++slot){
                if(oldIndices[startOld + slot] == neighbor){
                    alreadyPresent = true;
                    break;
                }
            }
            if(alreadyPresent){
                continue;
            }

            if(std::find(filtered.begin(), filtered.end(), neighbor) != filtered.end()){
                continue;
            }

            filtered.push_back(neighbor);
        }

        const int availableSlots = MAX_NEIGHBORS - oldCount;
        if(static_cast<int>(filtered.size()) > availableSlots){
            trimmedNeighborCount += static_cast<std::size_t>(filtered.size() - availableSlots);
            affectedAtomCount++;
            filtered.resize(static_cast<std::size_t>(availableSlots));
        }

        newOffsets[i + 1] = newOffsets[i] + oldCount + static_cast<int>(filtered.size());
    }

    if(trimmedNeighborCount > 0){
        spdlog::warn(
            "appendNeighbors trimmed {} appended neighbors across {} atoms to respect MAX_NEIGHBORS={}",
            trimmedNeighborCount,
            affectedAtomCount,
            static_cast<int>(MAX_NEIGHBORS)
        );
    }

    const size_t totalNeighbors = static_cast<size_t>(newOffsets[N]);
    auto newIndicesProp = std::make_shared<ParticleProperty>(
        totalNeighbors, DataType::Int, 1, 0, false);
    int* newIndices = newIndicesProp->dataInt();

    for(size_t i = 0; i < N; ++i){
        const int oldCount = _context.neighborCounts->getInt(static_cast<int>(i));
        const int startOld = oldOffsets[i];
        const int startNew = newOffsets[i];
        const auto& filtered = filteredExtras[i];
        std::copy(oldIndices + startOld, oldIndices + startOld + oldCount, newIndices + startNew);
        if(!filtered.empty()){
            std::copy(filtered.begin(), filtered.end(), newIndices + startNew + oldCount);
        }
        _context.neighborCounts->setInt(static_cast<int>(i), oldCount + static_cast<int>(filtered.size()));
    }

    _context.neighborOffsets = std::make_shared<ParticleProperty>(
        N + 1, DataType::Int, 1, 0, false);
    std::copy(newOffsets.begin(), newOffsets.end(), _context.neighborOffsets->dataInt());
    _context.neighborIndices = std::move(newIndicesProp);
}

void StructureAnalysis::setNeighborLatticeVectorOverrides(
    std::vector<Vector3> overrides,
    std::size_t stride
){
    _neighborLatticeVectorOverrides = std::move(overrides);
    _neighborLatticeVectorOverrideStride = stride;
}

void StructureAnalysis::setClusterRuleProvider(
    std::shared_ptr<const ClusterRuleProvider> clusterRuleProvider
){
    _clusterRuleProvider = std::move(clusterRuleProvider);
}

void StructureAnalysis::setCrystalInfoProvider(
    std::shared_ptr<const StructureAnalysisCrystalInfo> crystalInfoProvider
){
    _crystalInfoProvider = std::move(crystalInfoProvider);
}

AnalysisContext& StructureAnalysis::analysisContext(){
    return requireAnalysisContext(_context);
}

const AnalysisContext& StructureAnalysis::analysisContext() const{
    return requireAnalysisContext(_context);
}

int StructureAnalysis::findClosestSymmetryPermutation(int structureType, const Matrix3& rotation) const{
    return requireCrystalInfo(_crystalInfoProvider).findClosestSymmetryPermutation(structureType, rotation);
}

int StructureAnalysis::coordinationNumber(int structureType) const{
    return requireCrystalInfo(_crystalInfoProvider).coordinationNumber(structureType);
}

int StructureAnalysis::commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const{
    return requireCrystalInfo(_crystalInfoProvider).commonNeighborIndex(
        structureType,
        neighborIndex,
        commonNeighborSlot
    );
}

int StructureAnalysis::symmetryPermutationCount(int structureType) const{
    return requireCrystalInfo(_crystalInfoProvider).symmetryPermutationCount(structureType);
}

int StructureAnalysis::symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const{
    return requireCrystalInfo(_crystalInfoProvider).symmetryPermutationEntry(
        structureType,
        symmetryIndex,
        neighborIndex
    );
}

const Matrix3& StructureAnalysis::symmetryTransformation(int structureType, int symmetryIndex) const{
    return requireCrystalInfo(_crystalInfoProvider).symmetryTransformation(structureType, symmetryIndex);
}

int StructureAnalysis::symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const{
    return requireCrystalInfo(_crystalInfoProvider).symmetryInverseProduct(
        structureType,
        symmetryIndex,
        transformationIndex
    );
}

const Vector3& StructureAnalysis::latticeVector(int structureType, int latticeVectorIndex) const{
    return requireCrystalInfo(_crystalInfoProvider).latticeVector(structureType, latticeVectorIndex);
}

const Vector3& StructureAnalysis::neighborLatticeVector(int centralAtomIndex, int neighborIndex) const{
    if(hasNeighborLatticeVectorOverrides()){
        const std::size_t stride = _neighborLatticeVectorOverrideStride;
        assert(stride > 0);
        const std::size_t flatIndex =
            static_cast<std::size_t>(centralAtomIndex) * stride +
            static_cast<std::size_t>(neighborIndex);
        assert(flatIndex < _neighborLatticeVectorOverrides.size());
        return _neighborLatticeVectorOverrides[flatIndex];
    }

    const auto& context = analysisContext();
    assert(context.atomSymmetryPermutations);
    const int structureType = context.structureTypes->getInt(centralAtomIndex);
    assert(neighborIndex >= 0 && neighborIndex < coordinationNumber(structureType));
    const int symmetryPermutationIndex = context.atomSymmetryPermutations->getInt(centralAtomIndex);
    assert(symmetryPermutationIndex >= 0 && symmetryPermutationIndex < symmetryPermutationCount(structureType));
    return latticeVector(structureType, symmetryPermutationEntry(structureType, symmetryPermutationIndex, neighborIndex));
}

}
