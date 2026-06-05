#include <volt/core/reconstructed_structure.h>
#include <volt/analysis/reconstructed_dump_utils.h>

#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace Volt{

ReconstructedStructureContext::ReconstructedStructureContext(
    ParticleProperty* positions,
    const SimulationCell& cell
) :
    StructureContext(positions, cell){}

bool ReconstructedStructureContext::loadFromFrame(
    const LammpsParser::Frame& frame,
    StructureAnalysis& analysis,
    ReconstructedStructureContext& context,
    std::string* errorMessage
){
    if(frame.findAtomProperty("structure_type")){
        std::vector<int> structureTypes;
        if(!AnalysisDumpUtils::extractIntegralColumn(frame, "structure_type", structureTypes, errorMessage)){
            return false;
        }
        context._ownedStructureTypes = AnalysisDumpUtils::makeIntProperty(structureTypes);
        context.structureTypes = context._ownedStructureTypes.get();
    }

    std::vector<int> clusterIds;
    if(!AnalysisDumpUtils::extractIntegralColumn(frame, "cluster_id", clusterIds, errorMessage)){
        return false;
    }
    context.atomClusters = AnalysisDumpUtils::makeIntProperty(clusterIds);

    // Get direct pointers to neighbor index columns (avoid copying 18 × 7.3M ints)
    const size_t natoms = static_cast<size_t>(frame.natoms);
    std::array<const int*, MAX_NEIGHBORS> neighborSlotPtrs{};
    std::vector<std::vector<int>> neighborSlotStorage;

    for(int neighborSlot = 0; neighborSlot < MAX_NEIGHBORS; ++neighborSlot){
        const auto* column = frame.findAtomProperty(AnalysisDumpUtils::neighborIndexName(neighborSlot));
        if(!column){
            AnalysisDumpUtils::setError(errorMessage,
                "Missing required dump column '" + AnalysisDumpUtils::neighborIndexName(neighborSlot) + "'");
            return false;
        }
        if(column->dataType == DataType::Int && column->ints.size() == natoms){
            neighborSlotPtrs[neighborSlot] = column->ints.data();
        }else{
            std::vector<int> tmp;
            if(!AnalysisDumpUtils::extractIntegralColumn(frame,
                AnalysisDumpUtils::neighborIndexName(neighborSlot), tmp, errorMessage)){
                return false;
            }
            neighborSlotStorage.push_back(std::move(tmp));
            neighborSlotPtrs[neighborSlot] = neighborSlotStorage.back().data();
        }
    }

    // Parallel neighbor counting
    std::vector<int> neighborCounts(natoms, 0);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, natoms, 8192),
        [&](const tbb::blocked_range<size_t>& r){
        for(size_t atomIndex = r.begin(); atomIndex < r.end(); ++atomIndex){
            int count = 0;
            for(int s = 0; s < MAX_NEIGHBORS; ++s){
                if(neighborSlotPtrs[s][atomIndex] < 0) break;
                ++count;
            }
            neighborCounts[atomIndex] = count;
        }
    });

    // Prefix sum for offsets
    std::vector<int> compactOffsets(natoms + 1, 0);
    int totalNeighborEntries = 0;
    for(size_t i = 0; i < natoms; ++i){
        compactOffsets[i] = totalNeighborEntries;
        totalNeighborEntries += neighborCounts[i];
    }
    compactOffsets[natoms] = totalNeighborEntries;

    context.neighborCounts = AnalysisDumpUtils::makeIntProperty(neighborCounts);
    context.neighborOffsets = AnalysisDumpUtils::makeIntProperty(compactOffsets);

    // Parallel compact neighbor indices
    auto compactNeighborIndices = std::make_shared<ParticleProperty>(
        static_cast<size_t>(totalNeighborEntries), DataType::Int, 1, 0, true
    );
    int* compactData = compactNeighborIndices->dataInt();

    tbb::parallel_for(tbb::blocked_range<size_t>(0, natoms, 8192),
        [&](const tbb::blocked_range<size_t>& r){
        for(size_t atomIndex = r.begin(); atomIndex < r.end(); ++atomIndex){
            const int count = neighborCounts[atomIndex];
            const int start = compactOffsets[atomIndex];
            for(int s = 0; s < count; ++s){
                compactData[start + s] = neighborSlotPtrs[s][atomIndex];
            }
        }
    });

    context.neighborIndices = compactNeighborIndices;

    // Parallel max neighbor distance
    context.maximumNeighborDistance = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, natoms, 8192),
        0.0,
        [&](const tbb::blocked_range<size_t>& r, double localMax) -> double{
            for(size_t atomIndex = r.begin(); atomIndex < r.end(); ++atomIndex){
                const int count = neighborCounts[atomIndex];
                const int start = compactOffsets[atomIndex];
                for(int s = 0; s < count; ++s){
                    const int neighbor = compactData[start + s];
                    const Vector3 delta = context.simCell.wrapVector(
                        frame.positions[static_cast<size_t>(neighbor)] -
                        frame.positions[atomIndex]
                    );
                    localMax = std::max(localMax, delta.length());
                }
            }
            return localMax;
        },
        [](double a, double b){ return std::max(a, b); }
    );

    // Get direct pointers to lattice vector columns (avoid copying 54 × 7.3M doubles)
    const std::array<char, 3> axes = { 'x', 'y', 'z' };
    std::array<std::array<const double*, MAX_NEIGHBORS>, 3> latticePtrs{};
    std::vector<std::vector<double>> latticeStorage;

    for(int neighborSlot = 0; neighborSlot < MAX_NEIGHBORS; ++neighborSlot){
        for(size_t axisIndex = 0; axisIndex < 3; ++axisIndex){
            std::string colName = AnalysisDumpUtils::neighborLatticeComponentName(axes[axisIndex], neighborSlot);
            const auto* column = frame.findAtomProperty(colName);
            if(!column){
                AnalysisDumpUtils::setError(errorMessage, "Missing required dump column '" + colName + "'");
                return false;
            }
            if(column->dataType == DataType::Double && column->doubles.size() == natoms){
                latticePtrs[axisIndex][neighborSlot] = column->doubles.data();
            }else{
                std::vector<double> tmp;
                if(!AnalysisDumpUtils::extractDoubleColumn(frame, colName, tmp, errorMessage)){
                    return false;
                }
                latticeStorage.push_back(std::move(tmp));
                latticePtrs[axisIndex][neighborSlot] = latticeStorage.back().data();
            }
        }
    }

    // Parallel vector overrides assembly
    std::vector<Vector3> neighborVectorOverrides(
        context.atomCount() * static_cast<size_t>(MAX_NEIGHBORS),
        Vector3::Zero()
    );

    tbb::parallel_for(tbb::blocked_range<size_t>(0, context.atomCount(), 8192),
        [&](const tbb::blocked_range<size_t>& r){
        for(size_t atomIndex = r.begin(); atomIndex < r.end(); ++atomIndex){
            for(int s = 0; s < MAX_NEIGHBORS; ++s){
                neighborVectorOverrides[atomIndex * static_cast<size_t>(MAX_NEIGHBORS) + static_cast<size_t>(s)] = Vector3(
                    latticePtrs[0][s][atomIndex],
                    latticePtrs[1][s][atomIndex],
                    latticePtrs[2][s][atomIndex]
                );
            }
        }
    });

    analysis.setNeighborLatticeVectorOverrides(
        std::move(neighborVectorOverrides),
        static_cast<size_t>(MAX_NEIGHBORS)
    );

    return true;
}

bool ReconstructedStructureLoader::load(
    const LammpsParser::Frame& frame,
    const ClusterGraphExportPaths& paths,
    StructureAnalysis& structureAnalysis,
    ReconstructedStructureContext& context,
    std::string* errorMessage
){
    if(!ReconstructedStructureContext::loadFromFrame(frame, structureAnalysis, context, errorMessage)){
        return false;
    }

    if(!importClusterGraph(structureAnalysis, paths, errorMessage)){
        return false;
    }

    return true;
}

}
