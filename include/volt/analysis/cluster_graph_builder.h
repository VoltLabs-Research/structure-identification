#pragma once

#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/structure_analysis.h>

#include <deque>

namespace Volt{

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
    void formSuperClusters();
    void dissolveSmallClusters(int minClusterSize = kDefaultSmallClusterSize);
    void applyPreferredOrientation(Cluster* cluster);
    void reorientAtomsToAlign();
    void connectClusters();
    void build(bool dissolveSmallClusters = false);
    void grow(
        Cluster* cluster,
        std::deque<int>& atomsToVisit,
        Matrix_3<double>& orientationV,
        Matrix_3<double>& orientationW,
        int structureType
    );

protected:
    int selectInitialSymmetryPermutation(int atomIndex, int structureType) const;
    void buildClusterAssignments();

    AnalysisContext& _context;
    StructureAnalysis& _sa;
};

}
