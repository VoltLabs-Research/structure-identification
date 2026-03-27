#pragma once

#include <volt/analysis/analysis_context.h>
#include <volt/analysis/structure_analysis.h>

namespace Volt{

class ClusterBuilder{
public:
    static constexpr int kDefaultSmallClusterSize = 50;

    ClusterBuilder(
        StructureAnalysis& sa,
        AnalysisContext& context
    );

    Matrix3 quaternionToMatrix(const Quaternion& q) const;

    Cluster* startNew(int atomIndex, int structureType);
    Cluster* getParentGrain(Cluster* cluster);
    bool alreadyProcessedAtom(int index);
    void dissolveSmallClusters(int minClusterSize = kDefaultSmallClusterSize);

protected:
    AnalysisContext& _context;
    StructureAnalysis& _sa;
};

}
