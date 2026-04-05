#pragma once

#include <volt/analysis/structure_analysis.h>

namespace Volt::ClusterHierarchyUtils{

inline Cluster* getParentGrain(StructureAnalysis& structureAnalysis, Cluster* cluster){
    if(!cluster || !cluster->parentTransition){
        return cluster;
    }

    ClusterTransition* parentTransition = cluster->parentTransition;
    Cluster* parent = parentTransition->cluster2;
    while(parent && parent->parentTransition){
        parentTransition = structureAnalysis.clusterGraph().concatenateClusterTransitions(
            parentTransition,
            parent->parentTransition
        );
        parent = parent->parentTransition->cluster2;
    }

    cluster->parentTransition = parentTransition;
    return parent;
}

template<typename TransitionContainer>
void rebuildParentHierarchy(StructureAnalysis& structureAnalysis, const TransitionContainer& transitions){
    ClusterGraph& graph = structureAnalysis.clusterGraph();

    for(Cluster* cluster : graph.clusters()){
        if(!cluster || cluster->id == 0){
            continue;
        }
        cluster->parentTransition = nullptr;
        cluster->rank = 0;
    }

    for(ClusterTransition* transition : transitions){
        if(!transition || !transition->cluster1 || !transition->cluster2){
            continue;
        }

        Cluster* parent1 = getParentGrain(structureAnalysis, transition->cluster1);
        Cluster* parent2 = getParentGrain(structureAnalysis, transition->cluster2);
        if(!parent1 || !parent2 || parent1 == parent2){
            continue;
        }

        ClusterTransition* parentTransition = transition;
        if(parent2 != transition->cluster2){
            parentTransition = graph.concatenateClusterTransitions(
                parentTransition,
                transition->cluster2->parentTransition
            );
        }
        if(parent1 != transition->cluster1){
            parentTransition = graph.concatenateClusterTransitions(
                transition->cluster1->parentTransition->reverse,
                parentTransition
            );
        }

        if(parent1->rank > parent2->rank){
            parent2->parentTransition = parentTransition->reverse;
            continue;
        }

        parent1->parentTransition = parentTransition;
        if(parent1->rank == parent2->rank){
            parent2->rank++;
        }
    }

    for(Cluster* cluster : graph.clusters()){
        if(!cluster || cluster->id == 0){
            continue;
        }
        getParentGrain(structureAnalysis, cluster);
    }
}

}
