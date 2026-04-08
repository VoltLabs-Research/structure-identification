#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/structures/cluster_graph.h>

#include <string>

namespace Volt{

struct ClusterGraphExportPaths{
    std::string clustersTablePath;
    std::string clusterTransitionsTablePath;
};

void normalizeReconstructedClusterGraphForExport(
    StructureAnalysis& structureAnalysis,
    AnalysisContext& context
);

bool exportClusterGraph(
    StructureAnalysis& structureAnalysis,
    const AnalysisContext& context,
    const std::string& outputBase,
    ClusterGraphExportPaths* paths = nullptr
);

bool importClusterGraph(
    StructureAnalysis& structureAnalysis,
    const ClusterGraphExportPaths& paths,
    std::string* errorMessage = nullptr
);

void rebuildImportedClusterParentHierarchy(StructureAnalysis& structureAnalysis);

}
