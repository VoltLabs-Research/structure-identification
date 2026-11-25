#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/structures/cluster_graph.h>

#include <cstddef>
#include <iosfwd>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace Volt{

struct ClusterGraphExportPaths{
    std::string clustersTablePath;
    std::string clusterTransitionsTablePath;
};

// The rows exportClusterGraph() wrote to the clusters and cluster transitions tables.
struct ClusterGraphTables{
    std::size_t clusters = 0;
    std::size_t transitions = 0;
};

bool exportClusterGraph(
    ClusterGraph& clusterGraph,
    const std::string& outputBase,
    ClusterGraphExportPaths* paths = nullptr
);

bool importClusterGraph(
    StructureAnalysis& structureAnalysis,
    const ClusterGraphExportPaths& paths,
    std::string* errorMessage = nullptr
);

void rebuildImportedClusterParentHierarchy(StructureAnalysis& structureAnalysis);

// The table files and their rows, and the helpers of their export and import.
namespace ClusterGraphExportDetail{

inline constexpr const char* kClustersSuffix = "_clusters.table";
inline constexpr const char* kClusterTransitionsSuffix = "_cluster_transitions.table";
inline constexpr const char* kUnknownTopologyToken = "-";

struct ClusterRow{
    int clusterId = 0;
    std::string topologyName;
    Matrix3 orientation = Matrix3::Identity();
};

struct TransitionRow{
    int cluster1Id = 0;
    int cluster2Id = 0;
    Matrix3 tm = Matrix3::Identity();
    int distance = 1;
};

void writeMatrixColumns(std::ostream& stream, const Matrix3& matrix);
std::vector<Cluster*> sortedClusters(const ClusterGraph& clusterGraph);
std::vector<ClusterTransition*> collectExportedTransitions(const ClusterGraph& clusterGraph);
bool isTopologySymmetryEquivalent(std::string_view topologyName, const Matrix3& misorientation);
void rebuildParentHierarchy(StructureAnalysis& structureAnalysis);
std::vector<std::string> splitFields(const std::string& line);
bool readRows(
    const std::string& path,
    std::vector<std::string>& header,
    std::vector<std::vector<std::string>>& rows,
    std::string* errorMessage
);
bool requireColumnIndex(
    const std::unordered_map<std::string, std::size_t>& indices,
    const std::string& columnName,
    std::size_t& index,
    std::string* errorMessage
);
bool loadClustersTable(
    const std::string& path,
    std::vector<ClusterRow>& rows,
    std::string* errorMessage
);
bool loadTransitionsTable(
    const std::string& path,
    std::vector<TransitionRow>& rows,
    std::string* errorMessage
);
bool rebuildClusterGraph(
    StructureAnalysis& structureAnalysis,
    const std::vector<ClusterRow>& clusterRows,
    const std::vector<TransitionRow>& transitionRows,
    std::string* errorMessage
);
bool writeClustersTable(
    ClusterGraph& clusterGraph,
    const std::string& outputPath,
    std::size_t& rowCount
);
bool writeClusterTransitionsTable(
    const ClusterGraph& clusterGraph,
    const std::string& outputPath,
    std::size_t& rowCount
);

}

}
