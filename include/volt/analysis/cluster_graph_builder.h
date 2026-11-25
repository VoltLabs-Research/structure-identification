#pragma once

#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/structure_analysis.h>

#include <cstddef>
#include <deque>
#include <vector>

namespace Volt{

class ProgressReporter;

// What growing clusters from the seed atoms made, out of the crystalline atoms it started from. The atoms the
// clusters do not cover are seeds that had no allowed symmetry permutation, left unclustered as Other.
struct ClusterGrowth{
    std::size_t clusters = 0;
    std::size_t atoms = 0;
    std::size_t crystallineAtoms = 0;
};

// The cluster graph's transitions once the neighboring clusters are connected, and the neighbors connecting them
// left out so that no atom holds more than MAX_NEIGHBORS.
struct ClusterConnections{
    std::size_t transitions = 0;
    TrimmedNeighbors trimmedNeighbors;
};

// The transitions that join two clusters across a defect cluster between them.
struct SuperClusterBridges{
    std::size_t transitions = 0;
};

// The clusters marked as Other for having fewer atoms than minClusterSize.
struct DissolvedClusters{
    std::size_t clusters = 0;
    int minClusterSize = 0;
};

// The clusters build() left in the analysis's cluster graph, the cluster of each atom, and whether it dissolved the
// small ones. A caller that grew the clusters over some of the atoms only (the sites of a sublattice) names them, and
// the line then says how many of them are in a cluster.
struct BuiltClusters{
    const ClusterGraph& graph;
    const ParticleProperty& atomClusters;
    bool smallClustersDissolved = false;
    const std::vector<int>* latticeSites = nullptr;
};

class ClusterBuilder{
public:
    static constexpr int kDefaultSmallClusterSize = 50;

    ClusterBuilder(
        StructureAnalysis& sa,
        AnalysisContext& context
    );

    Cluster* startNew(int atomIndex, int structureType);
    Cluster* getParentGrain(Cluster* cluster);

    bool alreadyProcessedAtom(int index);
    bool calculateMisorientation(
        int atomIndex,
        int neighbor,
        int neighborIndex,
        Matrix3& outTransition
    );

    void processDefectCluster(Cluster* defectCluster);
    SuperClusterBridges formSuperClusters();
    DissolvedClusters dissolveSmallClusters(int minClusterSize = kDefaultSmallClusterSize);
    void applyPreferredOrientation(Cluster* cluster);
    void reorientAtomsToAlign();
    ClusterConnections connectClusters();
    BuiltClusters build(bool dissolveSmallClusters = false);
    void grow(
        Cluster* cluster,
        std::deque<int>& atomsToVisit,
        Matrix_3<double>& orientationV,
        Matrix_3<double>& orientationW,
        int structureType
    );

protected:
    int selectInitialSymmetryPermutation(int atomIndex, int structureType) const;
    ClusterGrowth buildClusterAssignments();

    AnalysisContext& _context;
    StructureAnalysis& _sa;

private:
    // grow() is public and its signature has no room for a reporter, so buildClusterAssignments()
    // lends it one for the duration of the flood fill.
    static thread_local ProgressReporter* _growthProgress;

    struct GrowthProgressScope{
        explicit GrowthProgressScope(ProgressReporter& progress);
        ~GrowthProgressScope();
    };
};

}
