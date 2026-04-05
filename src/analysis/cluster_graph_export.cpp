#include <volt/analysis/cluster_graph_export.h>
#include <volt/analysis/analysis_dump_utils.h>
#include <volt/analysis/cluster_hierarchy_utils.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Volt{

namespace ClusterGraphExportDetail{

constexpr const char* kClustersSuffix = "_clusters.table";
constexpr const char* kClusterTransitionsSuffix = "_cluster_transitions.table";

struct ClusterRow{
    int clusterId = 0;
    int structureType = 0;
    Matrix3 orientation = Matrix3::Identity();
};

struct TransitionRow{
    int cluster1Id = 0;
    int cluster2Id = 0;
    Matrix3 tm = Matrix3::Identity();
    int distance = 1;
};

void writeMatrixColumns(std::ostream& stream, const Matrix3& matrix){
    for(int row = 0; row < 3; ++row){
        for(int column = 0; column < 3; ++column){
            stream << '\t' << matrix(row, column);
        }
    }
}

std::vector<Cluster*> sortedClusters(const ClusterGraph& clusterGraph){
    std::vector<Cluster*> clusters;
    clusters.reserve(clusterGraph.clusters().size());

    for(Cluster* cluster : clusterGraph.clusters()){
        if(!cluster || cluster->id == 0){
            continue;
        }
        clusters.push_back(cluster);
    }

    std::sort(clusters.begin(), clusters.end(), [](const Cluster* left, const Cluster* right){
        return left->id < right->id;
    });
    return clusters;
}

std::vector<ClusterTransition*> collectDirectedTransitions(const ClusterGraph& clusterGraph){
    std::vector<ClusterTransition*> transitions;
    std::unordered_map<const ClusterTransition*, int> seen;

    for(Cluster* cluster : sortedClusters(clusterGraph)){
        for(ClusterTransition* transition = cluster->transitions; transition; transition = transition->next){
            if(!transition || !transition->cluster1 || !transition->cluster2){
                continue;
            }
            if(transition->cluster1->id == 0 || transition->cluster2->id == 0){
                continue;
            }
            if(seen.emplace(transition, static_cast<int>(seen.size())).second){
                transitions.push_back(transition);
            }
        }
    }

    return transitions;
}

void rebuildParentHierarchy(StructureAnalysis& structureAnalysis){
    ClusterGraph& graph = structureAnalysis.clusterGraph();
    std::vector<ClusterTransition*> superClusterTransitions;
    superClusterTransitions.reserve(graph.clusterTransitions().size());
    for(ClusterTransition* transition : graph.clusterTransitions()){
        if(transition && transition->distance == 2){
            superClusterTransitions.push_back(transition);
        }
    }
    ClusterHierarchyUtils::rebuildParentHierarchy(structureAnalysis, superClusterTransitions);
}

std::vector<std::string> splitFields(const std::string& line){
    std::istringstream stream(line);
    std::vector<std::string> fields;
    std::string field;
    while(stream >> field){
        fields.push_back(field);
    }
    return fields;
}

bool readRows(
    const std::string& path,
    std::vector<std::string>& header,
    std::vector<std::vector<std::string>>& rows,
    std::string* errorMessage
){
    std::ifstream input(path);
    if(!input.is_open()){
        AnalysisDumpUtils::setError(errorMessage, "Unable to open table file '" + path + "'");
        return false;
    }

    std::string line;
    while(std::getline(input, line)){
        header = splitFields(line);
        if(!header.empty()){
            break;
        }
    }
    if(header.empty()){
        AnalysisDumpUtils::setError(errorMessage, "Table file '" + path + "' is empty");
        return false;
    }

    rows.clear();
    while(std::getline(input, line)){
        auto fields = splitFields(line);
        if(fields.empty()){
            continue;
        }
        if(fields.size() != header.size()){
            AnalysisDumpUtils::setError(errorMessage, "Malformed row in table file '" + path + "'");
            return false;
        }
        rows.push_back(std::move(fields));
    }

    return true;
}

bool requireColumnIndex(
    const std::unordered_map<std::string, std::size_t>& indices,
    const std::string& columnName,
    std::size_t& index,
    std::string* errorMessage
){
    auto it = indices.find(columnName);
    if(it == indices.end()){
        AnalysisDumpUtils::setError(errorMessage, "Missing required table column '" + columnName + "'");
        return false;
    }
    index = it->second;
    return true;
}

bool loadClustersTable(
    const std::string& path,
    std::vector<ClusterRow>& rows,
    std::string* errorMessage
){
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rawRows;
    if(!readRows(path, header, rawRows, errorMessage)){
        return false;
    }

    std::unordered_map<std::string, std::size_t> indices;
    for(std::size_t i = 0; i < header.size(); ++i){
        indices.emplace(header[i], i);
    }

    std::size_t clusterIdIndex = 0;
    std::size_t structureTypeIndex = 0;
    if(!requireColumnIndex(indices, "cluster_id", clusterIdIndex, errorMessage) ||
       !requireColumnIndex(indices, "structure_type", structureTypeIndex, errorMessage)){
        return false;
    }

    std::array<std::size_t, 9> orientationIndices{};
    bool hasOrientationColumns = true;
    for(int row = 0; row < 3; ++row){
        for(int column = 0; column < 3; ++column){
            auto it = indices.find("orientation_" + std::to_string(row) + std::to_string(column));
            if(it == indices.end()){
                hasOrientationColumns = false;
                break;
            }
            orientationIndices[static_cast<std::size_t>(row * 3 + column)] = it->second;
        }
        if(!hasOrientationColumns){
            break;
        }
    }

    rows.clear();
    rows.reserve(rawRows.size());
    for(const auto& rawRow : rawRows){
        ClusterRow rowData;
        if(!AnalysisDumpUtils::tryParseInt(rawRow[clusterIdIndex], rowData.clusterId) ||
           !AnalysisDumpUtils::tryParseInt(rawRow[structureTypeIndex], rowData.structureType)){
            AnalysisDumpUtils::setError(errorMessage, "Failed to parse integral values in clusters table");
            return false;
        }

        if(hasOrientationColumns){
            for(int row = 0; row < 3; ++row){
                for(int column = 0; column < 3; ++column){
                    double value = 0.0;
                    if(!AnalysisDumpUtils::tryParseDouble(
                        rawRow[orientationIndices[static_cast<std::size_t>(row * 3 + column)]],
                        value
                    )){
                        AnalysisDumpUtils::setError(errorMessage, "Failed to parse cluster orientation matrix");
                        return false;
                    }
                    rowData.orientation(row, column) = value;
                }
            }
        }

        rows.push_back(rowData);
    }

    return true;
}

bool loadTransitionsTable(
    const std::string& path,
    std::vector<TransitionRow>& rows,
    std::string* errorMessage
){
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rawRows;
    if(!readRows(path, header, rawRows, errorMessage)){
        return false;
    }

    std::unordered_map<std::string, std::size_t> indices;
    for(std::size_t i = 0; i < header.size(); ++i){
        indices.emplace(header[i], i);
    }

    std::size_t cluster1Index = 0;
    std::size_t cluster2Index = 0;
    std::size_t distanceIndex = 0;
    if(!requireColumnIndex(indices, "cluster1_id", cluster1Index, errorMessage) ||
       !requireColumnIndex(indices, "cluster2_id", cluster2Index, errorMessage) ||
       !requireColumnIndex(indices, "distance", distanceIndex, errorMessage)){
        return false;
    }

    std::array<std::size_t, 9> tmIndices{};
    for(int row = 0; row < 3; ++row){
        for(int column = 0; column < 3; ++column){
            if(!requireColumnIndex(
                indices,
                "tm_" + std::to_string(row) + std::to_string(column),
                tmIndices[static_cast<std::size_t>(row * 3 + column)],
                errorMessage
            )){
                return false;
            }
        }
    }

    rows.clear();
    rows.reserve(rawRows.size());
    for(const auto& rawRow : rawRows){
        TransitionRow rowData;
        if(!AnalysisDumpUtils::tryParseInt(rawRow[cluster1Index], rowData.cluster1Id) ||
           !AnalysisDumpUtils::tryParseInt(rawRow[cluster2Index], rowData.cluster2Id) ||
           !AnalysisDumpUtils::tryParseInt(rawRow[distanceIndex], rowData.distance)){
            AnalysisDumpUtils::setError(errorMessage, "Failed to parse integral values in cluster transitions table");
            return false;
        }

        for(int row = 0; row < 3; ++row){
            for(int column = 0; column < 3; ++column){
                double value = 0.0;
                if(!AnalysisDumpUtils::tryParseDouble(
                    rawRow[tmIndices[static_cast<std::size_t>(row * 3 + column)]],
                    value
                )){
                    AnalysisDumpUtils::setError(errorMessage, "Failed to parse cluster transition matrix");
                    return false;
                }
                rowData.tm(row, column) = value;
            }
        }

        rows.push_back(rowData);
    }

    return true;
}

bool rebuildClusterGraph(
    StructureAnalysis& structureAnalysis,
    const std::vector<ClusterRow>& clusterRows,
    const std::vector<TransitionRow>& transitionRows,
    std::string* errorMessage
){
    ClusterGraph& graph = structureAnalysis.clusterGraph();

    std::vector<ClusterRow> sortedClusterRows = clusterRows;
    std::sort(sortedClusterRows.begin(), sortedClusterRows.end(), [](const ClusterRow& left, const ClusterRow& right){
        return left.clusterId < right.clusterId;
    });

    for(const ClusterRow& row : sortedClusterRows){
        if(row.clusterId <= 0){
            AnalysisDumpUtils::setError(errorMessage, "Cluster table contains invalid cluster_id");
            return false;
        }
        if(graph.findCluster(row.clusterId) != nullptr){
            AnalysisDumpUtils::setError(errorMessage, "Cluster table contains duplicate cluster_id");
            return false;
        }

        Cluster* cluster = graph.createCluster(row.structureType, row.clusterId);
        cluster->orientation = row.orientation;
        cluster->predecessor = nullptr;
        cluster->parentTransition = nullptr;
        cluster->rank = 0;
        cluster->transitions = nullptr;
    }

    std::vector<TransitionRow> sortedTransitionRows = transitionRows;
    std::sort(sortedTransitionRows.begin(), sortedTransitionRows.end(), [](const TransitionRow& left, const TransitionRow& right){
        if(left.distance != right.distance){
            return left.distance < right.distance;
        }
        if(left.cluster1Id != right.cluster1Id){
            return left.cluster1Id < right.cluster1Id;
        }
        return left.cluster2Id < right.cluster2Id;
    });

    for(const TransitionRow& row : sortedTransitionRows){
        Cluster* cluster1 = graph.findCluster(row.cluster1Id);
        Cluster* cluster2 = graph.findCluster(row.cluster2Id);
        if(!cluster1 || !cluster2){
            AnalysisDumpUtils::setError(errorMessage, "Cluster transitions reference unknown cluster ids");
            return false;
        }

        if(row.cluster1Id == row.cluster2Id || row.distance == 0){
            if(!row.tm.equals(Matrix3::Identity(), CA_TRANSITION_MATRIX_EPSILON)){
                AnalysisDumpUtils::setError(errorMessage, "Self cluster transitions must use the identity matrix");
                return false;
            }
            graph.createSelfTransition(cluster1);
            continue;
        }

        graph.createClusterTransition(cluster1, cluster2, row.tm, row.distance);
    }

    return true;
}

bool writeClustersTable(
    ClusterGraph& clusterGraph,
    const std::string& outputPath
){
    std::ofstream output(outputPath);
    if(!output.is_open()){
        return false;
    }

    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    output << "cluster_id\tstructure_type";
    for(int row = 0; row < 3; ++row){
        for(int column = 0; column < 3; ++column){
            output << "\torientation_" << row << column;
        }
    }
    output << '\n';

    for(Cluster* cluster : sortedClusters(clusterGraph)){
        output
            << cluster->id
            << '\t' << cluster->structure;
        writeMatrixColumns(output, cluster->orientation);
        output
            << '\n';
    }

    return output.good();
}

bool writeClusterTransitionsTable(
    const ClusterGraph& clusterGraph,
    const std::string& outputPath
){
    std::ofstream output(outputPath);
    if(!output.is_open()){
        return false;
    }

    const std::vector<ClusterTransition*> transitions = collectDirectedTransitions(clusterGraph);

    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    output << "cluster1_id\tcluster2_id";
    for(int row = 0; row < 3; ++row){
        for(int column = 0; column < 3; ++column){
            output << "\ttm_" << row << column;
        }
    }
    output << "\tdistance\n";

    for(std::size_t transitionIndex = 0; transitionIndex < transitions.size(); ++transitionIndex){
        const ClusterTransition* transition = transitions[transitionIndex];
        output
            << transition->cluster1->id
            << '\t' << transition->cluster2->id;
        writeMatrixColumns(output, transition->tm);
        output
            << '\t' << transition->distance
            << '\n';
    }

    return output.good();
}

}

using namespace ClusterGraphExportDetail;

bool exportClusterGraph(
    ClusterGraph& clusterGraph,
    const std::string& outputBase,
    ClusterGraphExportPaths* paths
){
    if(outputBase.empty()){
        if(paths){
            paths->clustersTablePath.clear();
            paths->clusterTransitionsTablePath.clear();
        }
        return true;
    }

    const std::string clustersTablePath = outputBase + kClustersSuffix;
    const std::string clusterTransitionsTablePath = outputBase + kClusterTransitionsSuffix;

    if(!writeClustersTable(clusterGraph, clustersTablePath)){
        return false;
    }
    if(!writeClusterTransitionsTable(clusterGraph, clusterTransitionsTablePath)){
        return false;
    }

    if(paths){
        paths->clustersTablePath = clustersTablePath;
        paths->clusterTransitionsTablePath = clusterTransitionsTablePath;
    }
    return true;
}

bool importClusterGraph(
    StructureAnalysis& structureAnalysis,
    const ClusterGraphExportPaths& paths,
    std::string* errorMessage
){
    if(paths.clustersTablePath.empty()){
        AnalysisDumpUtils::setError(errorMessage, "Missing --clusters-table input");
        return false;
    }
    if(paths.clusterTransitionsTablePath.empty()){
        AnalysisDumpUtils::setError(errorMessage, "Missing --clusters-transitions input");
        return false;
    }

    std::vector<ClusterRow> clusterRows;
    if(!loadClustersTable(paths.clustersTablePath, clusterRows, errorMessage)){
        return false;
    }

    std::vector<TransitionRow> transitionRows;
    if(!loadTransitionsTable(paths.clusterTransitionsTablePath, transitionRows, errorMessage)){
        return false;
    }

    return rebuildClusterGraph(structureAnalysis, clusterRows, transitionRows, errorMessage);
}

void rebuildImportedClusterParentHierarchy(StructureAnalysis& structureAnalysis){
    rebuildParentHierarchy(structureAnalysis);
}

}
