#pragma once

#include <volt/analysis/cluster_rule_provider.h>

#include <map>
#include <memory>
#include <vector>

namespace Volt {

struct OrientationClusterAtomState{
    Matrix3 orientation = Matrix3::Identity();
    bool valid = false;
    int preferredSymmetry = -1;
};

class OrientationClusterRuleProvider final : public ClusterRuleProvider{
public:
    explicit OrientationClusterRuleProvider(
        std::shared_ptr<const std::vector<OrientationClusterAtomState>> atomStates
    );

    void initializeClusterSeed(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        Cluster& cluster,
        int seedAtomIndex,
        int structureType
    ) const override;

    ClusterRuleDecision tryAssignNeighbor(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        const Cluster& cluster,
        int currentAtomIndex,
        int neighborAtomIndex,
        int neighborIndex,
        int structureType,
        int& outNeighborSymmetry
    ) const override;

    ClusterRuleDecision tryCalculateTransition(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        int atomIndex,
        int neighborAtomIndex,
        int neighborIndex,
        Matrix3& outTransition
    ) const override;

private:
    const OrientationClusterAtomState* stateFor(int atomIndex) const;

    static bool orthonormalizeMatrix(const Matrix3& input, Matrix3& output);
    static double matrixTrace(const Matrix3& matrix);
    static bool symmetryAllowed(std::uint64_t mask, int symmetryIndex);

    Vector3 transformedIdealVector(
        const StructureAnalysis& analysis,
        int structureType,
        int symmetryIndex,
        int slot
    ) const;

    int countMatchingLocalOverlap(
        const StructureAnalysis& analysis,
        int structureType,
        int currentSymmetry,
        int neighborSlot,
        int neighborSymmetry,
        int reverseSlot
    ) const;

    bool canonicalAlignmentTrace(
        const StructureAnalysis& analysis,
        int structureType,
        const OrientationClusterAtomState& currentState,
        int currentSymmetry,
        const OrientationClusterAtomState& neighborState,
        int neighborSymmetry,
        double& outTrace
    ) const;

    bool structureRequiresOrientationFallback(
        const StructureAnalysis& analysis,
        int structureType
    ) const;

    std::shared_ptr<const std::vector<OrientationClusterAtomState>> _atomStates;
    mutable std::map<int, bool> _fallbackStructures;
};

}
