#include <volt/core/volt.h>
#include <volt/analysis/structure_analysis.h>

namespace Volt {

StructureAnalysis::StructureAnalysis(
    AnalysisContext& context,
    bool identifyPlanarDefects, 
    Mode identificationMode,
    float rmsd
) :
    _context(context),
    _identificationMode(identificationMode),
    _identifyPlanarDefects(identifyPlanarDefects),
    _rmsd(rmsd),    
    _clusterGraph(std::make_unique<ClusterGraph>())
{
    _context.atomSymmetryPermutations = std::make_shared<ParticleProperty>(
        _context.atomCount(), DataType::Int, 1, 0, true);

    _context.neighborOffsets = std::make_shared<ParticleProperty>(
        _context.atomCount() + 1, DataType::Int, 1, 0, true);

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

StructureAnalysis::~StructureAnalysis() = default;

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

}
