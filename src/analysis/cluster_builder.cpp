#include <volt/analysis/cluster_builder.h>

namespace Volt{

ClusterBuilder::ClusterBuilder(
    StructureAnalysis& sa,
    AnalysisContext& context
): _sa(sa), _context(context){}

Matrix3 ClusterBuilder::quaternionToMatrix(const Quaternion& q) const{
    double w = q.w(), x = q.x(), y = q.y(), z = q.z();
    Matrix3 R;
    R(0,0) = 1 - 2 * (y * y + z * z);
    R(0,1) = 2 * (x * y - w * z);
    R(0,2) = 2 * (x * z + w * y);
    R(1,0) = 2 * (x * y + w * z);
    R(1,1) = 1 - 2 * (x * x + z * z);
    R(1,2) = 2 * (y * z - w * x);
    R(2,0) = 2 * (x * z - w * y);
    R(2,1) = 2 * (y * z + w * x);
    R(2,2) = 1 - 2 * (x * x + y * y);
    return R;
}

bool ClusterBuilder::alreadyProcessedAtom(int index){
    return _context.atomClusters->getInt(index) != 0
            || _context.structureTypes->getInt(index) == StructureType::OTHER;
}

void ClusterBuilder::dissolveSmallClusters(int minClusterSize){
    for(Cluster* cluster : _sa.clusterGraph().clusters()){
        if(!cluster){
            continue;
        }

        int clusterSize = cluster->atomCount;
        if(cluster->parentTransition != nullptr && cluster->parentTransition->cluster2 != nullptr){
            clusterSize = cluster->parentTransition->cluster2->atomCount;
        }

        if(cluster->structure == LATTICE_OTHER || clusterSize >= minClusterSize){
            continue;
        }

        cluster->structure = LATTICE_OTHER;
    }
}

Cluster* ClusterBuilder::startNew(int atomIndex, int structureType){
    Cluster* cluster = _sa.clusterGraph().createCluster(structureType);
    assert(cluster->id > 0);

    cluster->atomCount = 1;
    _context.atomClusters->setInt(atomIndex, cluster->id);
    return cluster;
}

Cluster* ClusterBuilder::getParentGrain(Cluster* cluster){
    if(!cluster->parentTransition){
        return cluster;
    }

    ClusterTransition* parentTransition = cluster->parentTransition;
    Cluster* parent = parentTransition->cluster2;

    while(parent->parentTransition){
        parentTransition = _sa.clusterGraph().concatenateClusterTransitions(parentTransition, parent->parentTransition);
        parent = parent->parentTransition->cluster2;
    }

    cluster->parentTransition = parentTransition;
    return parent;
}

}
