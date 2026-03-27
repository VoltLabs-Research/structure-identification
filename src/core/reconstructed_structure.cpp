#include <volt/core/reconstructed_structure.h>
#include <volt/analysis/analysis_dump_utils.h>

#include <algorithm>
#include <array>
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
    std::string text;
    if(!AnalysisDumpUtils::extractHeaderValue(
        frame,
        AnalysisDumpUtils::kMaximumNeighborDistanceHeader,
        text,
        errorMessage
    )){
        return false;
    }
    if(!AnalysisDumpUtils::tryParseDouble(text, context.maximumNeighborDistance)){
        AnalysisDumpUtils::setError(
            errorMessage,
            "Failed to parse dump header '" + std::string(AnalysisDumpUtils::kMaximumNeighborDistanceHeader) + "'"
        );
        return false;
    }

    std::vector<int> clusterIds;
    if(!AnalysisDumpUtils::extractIntegralColumn(frame, "cluster_id", clusterIds, errorMessage)){
        return false;
    }
    context.atomClusters = AnalysisDumpUtils::makeIntProperty(clusterIds);

    std::vector<int> neighborCounts;
    if(!AnalysisDumpUtils::extractIntegralColumn(frame, "neighbor_counts", neighborCounts, errorMessage)){
        return false;
    }
    context.neighborCounts = AnalysisDumpUtils::makeIntProperty(neighborCounts);

    std::vector<std::vector<int>> expandedNeighborSlots(static_cast<std::size_t>(MAX_NEIGHBORS));
    for(int neighborSlot = 0; neighborSlot < MAX_NEIGHBORS; ++neighborSlot){
        if(!AnalysisDumpUtils::extractIntegralColumn(
            frame,
            AnalysisDumpUtils::neighborIndexName(neighborSlot),
            expandedNeighborSlots[static_cast<std::size_t>(neighborSlot)],
            errorMessage
        )){
            return false;
        }
    }

    std::vector<int> compactOffsets(static_cast<std::size_t>(frame.natoms) + 1, 0);
    int totalNeighborEntries = 0;
    for(int atomIndex = 0; atomIndex < frame.natoms; ++atomIndex){
        const int count = neighborCounts[static_cast<std::size_t>(atomIndex)];
        if(count < 0 || count > MAX_NEIGHBORS){
            AnalysisDumpUtils::setError(
                errorMessage,
                "neighbor_counts contains invalid value " + std::to_string(count) +
                " for atom index " + std::to_string(atomIndex)
            );
            return false;
        }

        compactOffsets[static_cast<std::size_t>(atomIndex)] = totalNeighborEntries;
        totalNeighborEntries += count;
    }
    compactOffsets[static_cast<std::size_t>(frame.natoms)] = totalNeighborEntries;

    context.neighborOffsets = AnalysisDumpUtils::makeIntProperty(compactOffsets);

    auto compactNeighborIndices = std::make_shared<ParticleProperty>(
        static_cast<std::size_t>(totalNeighborEntries), DataType::Int, 1, 0, true
    );
    int* compactData = compactNeighborIndices->dataInt();
    for(int atomIndex = 0; atomIndex < frame.natoms; ++atomIndex){
        const int count = neighborCounts[static_cast<std::size_t>(atomIndex)];
        const int start = compactOffsets[static_cast<std::size_t>(atomIndex)];
        for(int neighborSlot = 0; neighborSlot < count; ++neighborSlot){
            const int neighbor = expandedNeighborSlots[static_cast<std::size_t>(neighborSlot)][static_cast<std::size_t>(atomIndex)];
            if(neighbor < 0 || neighbor >= frame.natoms){
                AnalysisDumpUtils::setError(
                    errorMessage,
                    AnalysisDumpUtils::neighborIndexName(neighborSlot) +
                    " contains invalid atom index " + std::to_string(neighbor)
                );
                return false;
            }
            compactData[start + neighborSlot] = neighbor;
        }
    }
    context.neighborIndices = compactNeighborIndices;

    const std::array<char, 3> axes = { 'x', 'y', 'z' };
    std::array<std::vector<std::vector<double>>, 3> neighborLattice;
    for(auto& axisValues : neighborLattice){
        axisValues.resize(static_cast<std::size_t>(MAX_NEIGHBORS));
    }
    for(int neighborSlot = 0; neighborSlot < MAX_NEIGHBORS; ++neighborSlot){
        for(std::size_t axisIndex = 0; axisIndex < axes.size(); ++axisIndex){
            if(!AnalysisDumpUtils::extractDoubleColumn(
                frame,
                AnalysisDumpUtils::neighborLatticeComponentName(axes[axisIndex], neighborSlot),
                neighborLattice[axisIndex][static_cast<std::size_t>(neighborSlot)],
                errorMessage
            )){
                return false;
            }
        }
    }

    std::vector<Vector3> neighborVectorOverrides(
        context.atomCount() * static_cast<std::size_t>(MAX_NEIGHBORS),
        Vector3::Zero()
    );
    for(std::size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
        for(int neighborSlot = 0; neighborSlot < MAX_NEIGHBORS; ++neighborSlot){
            neighborVectorOverrides[atomIndex * static_cast<std::size_t>(MAX_NEIGHBORS) +
                                    static_cast<std::size_t>(neighborSlot)] = Vector3(
                neighborLattice[0][static_cast<std::size_t>(neighborSlot)][atomIndex],
                neighborLattice[1][static_cast<std::size_t>(neighborSlot)][atomIndex],
                neighborLattice[2][static_cast<std::size_t>(neighborSlot)][atomIndex]
            );
        }
    }
    analysis.setNeighborLatticeVectorOverrides(
        std::move(neighborVectorOverrides),
        static_cast<std::size_t>(MAX_NEIGHBORS)
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
