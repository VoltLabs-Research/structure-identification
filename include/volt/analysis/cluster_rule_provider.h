#pragma once

#include <volt/core/volt.h>

#include <memory>

namespace Volt{

class AnalysisContext;
class Cluster;
class StructureAnalysis;

enum class ClusterRuleDecision{
    Unhandled,
    Rejected,
    Accepted
};

class ClusterRuleProvider{
public:
    virtual ~ClusterRuleProvider() = default;

    virtual void initializeClusterSeed(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        Cluster& cluster,
        int seedAtomIndex,
        int structureType
    ) const{}

    virtual bool finalizeClusterOrientation(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        Cluster& cluster,
        int seedAtomIndex,
        int structureType
    ) const{
        return false;
    }

    virtual ClusterRuleDecision tryAssignNeighbor(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        const Cluster& cluster,
        int currentAtomIndex,
        int neighborAtomIndex,
        int neighborIndex,
        int structureType,
        int& outNeighborSymmetry
    ) const{
        return ClusterRuleDecision::Unhandled;
    }

    virtual ClusterRuleDecision tryCalculateTransition(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        int atomIndex,
        int neighborAtomIndex,
        int neighborIndex,
        Matrix3& outTransition
    ) const{
        return ClusterRuleDecision::Unhandled;
    }
};

}
