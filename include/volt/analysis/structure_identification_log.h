#pragma once

#include <volt/core/run_log.h>

namespace Volt{

class ReconstructedStructureContext;
struct AnnotatedDump;
struct BuiltClusters;
struct ClusterConnections;
struct ClusterGraphImport;
struct ClusterGraphTables;
struct ClusterGrowth;
struct DissolvedClusters;
struct NeighborLatticeTable;
struct NeighborShellOrder;
struct SuperClusterBridges;

namespace AnalysisPipelineUtils{
struct ClusterOutputs;
}

// The execution log of a structure identification. The library writes the indented lines of its own passes
// through one, and the log of each plugin built on the library derives from it for the lines they all print: the
// clusters and the cluster outputs. Scripts parse these lines, so their wording stays.
class StructureIdentificationLog : public RunLog{
public:
    using RunLog::dump;

    void dump(const ClusterGrowth& growth);
    void dump(const ClusterConnections& connections);
    void dump(const SuperClusterBridges& bridges);
    void dump(const DissolvedClusters& dissolved);
    void dump(const NeighborShellOrder& order);
    void dump(const AnnotatedDump& annotatedDump);
    void dump(const ClusterGraphTables& tables);
    void dump(const NeighborLatticeTable& neighborLattice);
    void dump(const ReconstructedStructureContext& context);
    void dump(const ClusterGraphImport& graphImport);

    void dump(const BuiltClusters& built);
    void dump(const AnalysisPipelineUtils::ClusterOutputs& outputs);
};

}
