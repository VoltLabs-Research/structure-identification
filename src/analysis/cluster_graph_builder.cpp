#include <volt/analysis/cluster_graph_builder.h>
#include <volt/analysis/cluster_hierarchy_rebuilder.h>
#include <volt/analysis/cluster_rule_provider.h>

#include <spdlog/spdlog.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include <array>
#include <cstdint>
#include <limits>

namespace Volt{

namespace ClusterBuilderDetail{

constexpr bool symmetryAllowed(std::uint64_t mask, int symmetryIndex){
    return symmetryIndex >= 0
        && symmetryIndex < 63
        && (mask & (std::uint64_t{1} << symmetryIndex)) != 0;
}

}

using namespace ClusterBuilderDetail;

ClusterBuilder::ClusterBuilder(
    StructureAnalysis& sa,
    AnalysisContext& context
): _sa(sa), _context(context){}

bool ClusterBuilder::alreadyProcessedAtom(int index){
    return _context.atomClusters->getInt(index) != 0
            || _context.structureTypes->getInt(index) == StructureType::OTHER;
}

int ClusterBuilder::selectInitialSymmetryPermutation(int atomIndex, int structureType) const{
    const std::uint64_t allowedMask = static_cast<std::uint64_t>(
        _context.atomAllowedSymmetryMasks->getInt64(static_cast<std::size_t>(atomIndex))
    );

    const int symmetryCount = _sa.symmetryPermutationCount(structureType);
    const int existingSymmetry = _context.atomSymmetryPermutations->getInt(static_cast<std::size_t>(atomIndex));
    if(existingSymmetry >= 0 &&
       existingSymmetry < symmetryCount &&
       symmetryAllowed(allowedMask, existingSymmetry)){
        return existingSymmetry;
    }

    for(int symmetryIndex = 0; symmetryIndex < symmetryCount; ++symmetryIndex){
        if(symmetryAllowed(allowedMask, symmetryIndex)){
            return symmetryIndex;
        }
    }

    return -1;
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
        cluster->topologyName.clear();
    }
}

void ClusterBuilder::grow(
    Cluster* cluster,
    std::deque<int>& atomsToVisit,
    Matrix_3<double>& orientationV,
    Matrix_3<double>& orientationW,
    int structureType
){
    const int coordinationNumber = _sa.coordinationNumber(structureType);

    while(!atomsToVisit.empty()){
        int currentAtomIndex = atomsToVisit.front();
        atomsToVisit.pop_front();

        const int currentSymmetryPermutation = _context.atomSymmetryPermutations->getInt(static_cast<size_t>(currentAtomIndex));
        if(currentSymmetryPermutation < 0){
            continue;
        }
        
        for(int neighborIndex = 0; neighborIndex < coordinationNumber; neighborIndex++){
            const int neighborAtomIndex = _sa.getNeighbor(currentAtomIndex, neighborIndex);
            if(neighborAtomIndex < 0 || neighborAtomIndex == currentAtomIndex){
                continue;
            }

            const int latticeVectorIndex = _sa.symmetryPermutationEntry(structureType, currentSymmetryPermutation, neighborIndex);
            const Vector3& latticeVector = _sa.latticeVector(structureType, latticeVectorIndex);
            const Vector3& spatialVector = _context.simCell.wrapVector(
                _context.positions->getPoint3(static_cast<size_t>(neighborAtomIndex)) -
                _context.positions->getPoint3(static_cast<size_t>(currentAtomIndex))
            );

            for(size_t i = 0; i < 3; i++){
                for(size_t j = 0; j < 3; j++){
                    orientationV(i, j) += latticeVector[j] * latticeVector[i];
                    orientationW(i, j) += latticeVector[j] * spatialVector[i];
                }
            }

            if(_context.atomClusters->getInt(neighborAtomIndex) != 0){
                continue;
            }

            const int neighborStructureType = _context.structureTypes->getInt(neighborAtomIndex);
            if(neighborStructureType != structureType){
                const auto* crossProvider = _sa.clusterRuleProvider();
                if(neighborStructureType == LATTICE_OTHER ||
                   !crossProvider ||
                   !crossProvider->allowsCrossStructureGrowth()){
                    continue;
                }
            }

            if(const auto* clusterRuleProvider = _sa.clusterRuleProvider()){
                int resolvedNeighborSymmetry = -1;
                switch(clusterRuleProvider->tryAssignNeighbor(
                    _sa,
                    _context,
                    *cluster,
                    currentAtomIndex,
                    neighborAtomIndex,
                    neighborIndex,
                    structureType,
                    resolvedNeighborSymmetry
                )){
                    case ClusterRuleDecision::Accepted:
                        _context.atomClusters->setInt(neighborAtomIndex, cluster->id);
                        cluster->atomCount++;
                        _context.atomSymmetryPermutations->setInt(neighborAtomIndex, resolvedNeighborSymmetry);
                        atomsToVisit.push_back(neighborAtomIndex);
                        continue;
                    case ClusterRuleDecision::Rejected:
                        continue;
                    case ClusterRuleDecision::Unhandled:
                        break;
                }
            }

            std::array<int, 3> matchedAtomIndices = { -1, -1, -1 };
            Matrix3 currentLocalBasis;
            bool properOverlap = true;

            for(int axis = 0; axis < 3; axis++){
                if(axis != 2){
                    int commonNeighbor = _sa.commonNeighborIndex(structureType, neighborIndex, axis);
                    if(commonNeighbor < 0){
                        properOverlap = false;
                        break;
                    }

                    matchedAtomIndices[static_cast<std::size_t>(axis)] = _sa.getNeighbor(currentAtomIndex, commonNeighbor);
                    if(matchedAtomIndices[static_cast<std::size_t>(axis)] < 0){
                        properOverlap = false;
                        break;
                    }
                    
                    const int commonNeighborSlot = _sa.symmetryPermutationEntry(
                        structureType,
                        currentSymmetryPermutation,
                        commonNeighbor
                    );

                    const int bondNeighborSlot = _sa.symmetryPermutationEntry(
                        structureType,
                        currentSymmetryPermutation,
                        neighborIndex
                    );

                    const Vector3 commonNeighborVector = _sa.latticeVector(structureType, commonNeighborSlot);
                    const Vector3 bondNeighborVector = _sa.latticeVector(structureType, bondNeighborSlot);
                    currentLocalBasis.column(axis) = commonNeighborVector - bondNeighborVector;
                }else{
                    matchedAtomIndices[static_cast<std::size_t>(axis)] = currentAtomIndex;
                    const int bondNeighborSlot = _sa.symmetryPermutationEntry(
                        structureType,
                        currentSymmetryPermutation,
                        neighborIndex
                    );
                    const Vector3 bondVector = _sa.latticeVector(structureType, bondNeighborSlot);
                    currentLocalBasis.column(axis) = -bondVector;
                }
            }

            if(!properOverlap) continue;

            const std::uint64_t allowedMask = static_cast<std::uint64_t>(
                _context.atomAllowedSymmetryMasks->getInt64(static_cast<std::size_t>(neighborAtomIndex))
            );

            Matrix3 neighborCanonicalBasis;
            bool properReverseOverlap = true;
            for(int axis = 0; axis < 3; ++axis){
                const int reverseSlot = _sa.findNeighbor(
                    neighborAtomIndex,
                    matchedAtomIndices[static_cast<std::size_t>(axis)]
                );
                if(reverseSlot < 0){
                    properReverseOverlap = false;
                    break;
                }

                neighborCanonicalBasis.column(axis) = _sa.latticeVector(structureType, reverseSlot);
            }

            if(!properReverseOverlap){
                continue;
            }

            Matrix3 neighborCanonicalBasisInverse;
            if(!neighborCanonicalBasis.inverse(neighborCanonicalBasisInverse)){
                continue;
            }

            const Matrix3 transition = currentLocalBasis * neighborCanonicalBasisInverse;
            const int neighborSymmetryCount = _sa.symmetryPermutationCount(structureType);
            for(int symmetryIndex = 0; symmetryIndex < neighborSymmetryCount; ++symmetryIndex){
                if(!symmetryAllowed(allowedMask, symmetryIndex)){
                    continue;
                }
                if(!transition.equals(_sa.symmetryTransformation(structureType, symmetryIndex), CA_TRANSITION_MATRIX_EPSILON)){
                    continue;
                }

                _context.atomClusters->setInt(neighborAtomIndex, cluster->id);
                cluster->atomCount++;
                _context.atomSymmetryPermutations->setInt(neighborAtomIndex, symmetryIndex);
                atomsToVisit.push_back(neighborAtomIndex);
                break;
            }
        }
    }
}

void ClusterBuilder::applyPreferredOrientation(Cluster* cluster){
    double smallestDeviation = std::numeric_limits<double>::max();
    Matrix3 oldOrientation = cluster->orientation;

    for(int symmetryIndex = 0; symmetryIndex < _sa.symmetryPermutationCount(cluster->structure); symmetryIndex++){
        Matrix3 newOrientation = oldOrientation * _sa.symmetryTransformation(cluster->structure, symmetryIndex).inverse();
        double scaling = std::pow(std::abs(newOrientation.determinant()), 1.0 / 3.0);

        for(const auto& preferredOrientation : _context.preferredCrystalOrientations){
            double deviation = 0;
            for(size_t i = 0; i < 3; i++){
                for(size_t j = 0; j < 3; j++){
                    deviation += std::abs(newOrientation(i, j) / scaling - preferredOrientation(i, j));
                }
            }

            if(deviation < smallestDeviation){
                smallestDeviation = deviation;
                cluster->symmetryTransformation = symmetryIndex;
                cluster->orientation = newOrientation;
            }
        }
    }
}

bool ClusterBuilder::calculateMisorientation(
    int atomIndex,
    int neighbor,
    int neighborIndex,
    Matrix3& outTransition
){
    if(const auto* clusterRuleProvider = _sa.clusterRuleProvider()){
        switch(clusterRuleProvider->tryCalculateTransition(
            _sa,
            _context,
            atomIndex,
            neighbor,
            neighborIndex,
            outTransition
        )){
            case ClusterRuleDecision::Accepted:
                return true;
            case ClusterRuleDecision::Rejected:
                return false;
            case ClusterRuleDecision::Unhandled:
                break;
        }
    }

    int structureType = _context.structureTypes->getInt(atomIndex);
    const int neighborStructureType = _context.structureTypes->getInt(neighbor);
    int symIndex = _context.atomSymmetryPermutations->getInt(atomIndex);
    const int neighborSymIndex = _context.atomSymmetryPermutations->getInt(neighbor);
    if(symIndex < 0 || neighborSymIndex < 0){
        return false;
    }

    Matrix3 tm1, tm2;
    for(int i = 0; i < 3; i++){
        int ai;
        if(i != 2){
            int cnIdx = _sa.commonNeighborIndex(structureType, neighborIndex, i);
            if(cnIdx < 0){
                return false;
            }

            ai = _sa.getNeighbor(atomIndex, cnIdx);
            if(ai < 0){
                return false;
            }
            tm1.column(i) = _sa.latticeVector(structureType, _sa.symmetryPermutationEntry(structureType, symIndex, cnIdx)) -
                            _sa.latticeVector(structureType, _sa.symmetryPermutationEntry(structureType, symIndex, neighborIndex));
        }else{
            ai = atomIndex;
            tm1.column(i) = -_sa.latticeVector(structureType, _sa.symmetryPermutationEntry(structureType, symIndex, neighborIndex));
        }

        if(_sa.numberOfNeighbors(neighbor) != _sa.coordinationNumber(structureType)){
            return false;
        }

        int j = _sa.findNeighbor(neighbor, ai);
        if(j == -1){
            return false;
        }

        tm2.column(i) = _sa.latticeVector(
            neighborStructureType,
            _sa.symmetryPermutationEntry(neighborStructureType, neighborSymIndex, j)
        );
    }

    if(std::abs(tm1.determinant()) < EPSILON){
        return false;
    }

    Matrix3 tm1inv;
    if(!tm1.inverse(tm1inv)){
        return false;
    }

    outTransition = tm2 * tm1inv;
    return true;
}

void ClusterBuilder::connectClusters(){
    std::vector<std::vector<int>> extras(_context.atomCount());

    for(size_t atomIndex = 0; atomIndex < _context.atomCount(); ++atomIndex){
        int clusterId = _context.atomClusters->getInt(atomIndex);
        if(clusterId == 0){
            continue;
        }

        Cluster* cluster1 = _sa.clusterGraph().findCluster(clusterId);
        const int nn = _sa.numberOfNeighbors(atomIndex);

        for(int ni = 0; ni < nn; ++ni){
            int neighbor = _sa.getNeighbor(atomIndex, ni);
            if(neighbor < 0 || neighbor == static_cast<int>(atomIndex)){
                continue;
            }

            int neighborClusterId = _context.atomClusters->getInt(neighbor);
            if(neighborClusterId == 0){
                extras[neighbor].push_back(static_cast<int>(atomIndex));
                continue;
            }

            if(neighborClusterId == cluster1->id){
                continue;
            }

            Cluster* cluster2 = _sa.clusterGraph().findCluster(neighborClusterId);
            if(ClusterTransition* existing = cluster1->findTransition(cluster2)){
                existing->area++;
                existing->reverse->area++;
                continue;
            }

            Matrix3 transition;
            if(!calculateMisorientation(static_cast<int>(atomIndex), neighbor, ni, transition)){
                continue;
            }
            if(!transition.isOrthogonalMatrix()){
                continue;
            }

            ClusterTransition* transitionLink = _sa.clusterGraph().createClusterTransition(cluster1, cluster2, transition);
            transitionLink->area++;
            transitionLink->reverse->area++;
        }
    }

    _sa.appendNeighbors(extras);
    spdlog::info("Number of cluster transitions: {}", _sa.clusterGraph().clusterTransitions().size());
}

void ClusterBuilder::reorientAtomsToAlign(){
    tbb::parallel_for(tbb::blocked_range<size_t>(0, _context.atomCount()), [this](const tbb::blocked_range<size_t>& r){
        for(size_t atomIndex = r.begin(); atomIndex != r.end(); ++atomIndex){
            int clusterId = _context.atomClusters->getInt(atomIndex);
            if(clusterId == 0){
                continue;
            }

            Cluster* cluster = _sa.clusterGraph().findCluster(clusterId);
            assert(cluster);
            if(cluster->symmetryTransformation == 0){
                continue;
            }

            int oldSymmetry = _context.atomSymmetryPermutations->getInt(atomIndex);
            if(oldSymmetry < 0){
                continue;
            }
            int newSymmetry = _sa.symmetryInverseProduct(
                cluster->structure,
                oldSymmetry,
                cluster->symmetryTransformation
            );

            _context.atomSymmetryPermutations->setInt(atomIndex, newSymmetry);
        }
    });
}

void ClusterBuilder::buildClusterAssignments(){
    for(size_t seedAtomIndex = 0; seedAtomIndex < _context.atomCount(); seedAtomIndex++){
        if(alreadyProcessedAtom(seedAtomIndex)) continue;

        int structureType = _context.structureTypes->getInt(seedAtomIndex);
        const int seedSymmetryPermutation = selectInitialSymmetryPermutation(static_cast<int>(seedAtomIndex), structureType);
        if(seedSymmetryPermutation < 0){
            _context.structureTypes->setInt(seedAtomIndex, LATTICE_OTHER);
            continue;
        }

        Cluster* cluster = startNew(seedAtomIndex, structureType);
        _context.atomSymmetryPermutations->setInt(seedAtomIndex, seedSymmetryPermutation);
        if(const auto* clusterRuleProvider = _sa.clusterRuleProvider()){
            clusterRuleProvider->initializeClusterSeed(
                _sa,
                _context,
                *cluster,
                static_cast<int>(seedAtomIndex),
                structureType
            );
        }

        Matrix_3<double> orientationV = Matrix_3<double>::Zero();
        Matrix_3<double> orientationW = Matrix_3<double>::Zero();
        std::deque<int> atomsToVisit(1, seedAtomIndex);

        grow(cluster, atomsToVisit, orientationV, orientationW, structureType);

        bool clusterOrientationHandled = false;
        if(const auto* clusterRuleProvider = _sa.clusterRuleProvider()){
            clusterOrientationHandled = clusterRuleProvider->finalizeClusterOrientation(
                _sa,
                _context,
                *cluster,
                static_cast<int>(seedAtomIndex),
                structureType
            );
        }

        if(!clusterOrientationHandled){
            Matrix3 orientationVInverse;
            if(orientationV.inverse(orientationVInverse)){
                cluster->orientation = Matrix3(orientationW * orientationVInverse);
            }else{
                cluster->orientation = Matrix3::Identity();
            }

            if(structureType == _context.inputCrystalType && !_context.preferredCrystalOrientations.empty()){
                applyPreferredOrientation(cluster);
            }
        }
    }
}

void ClusterBuilder::build(bool dissolveSmallClusters){
    buildClusterAssignments();
    reorientAtomsToAlign();
    connectClusters();
    formSuperClusters();
    if(dissolveSmallClusters){
        ClusterBuilder::dissolveSmallClusters();
    }
}

void ClusterBuilder::processDefectCluster(Cluster* defectCluster){
    for(ClusterTransition* transition = defectCluster->transitions; transition; transition = transition->next){
        if(transition->cluster2->structure != _context.inputCrystalType || transition->distance != 1){
            continue;
        }
        for(ClusterTransition* sibling = transition->next; sibling; sibling = sibling->next) {
            if(sibling->cluster2->structure != _context.inputCrystalType || sibling->distance != 1){
                continue;
            }
            if(sibling->cluster2 == transition->cluster2){
                continue;
            }

            Matrix3 misorientation = sibling->tm * transition->reverse->tm;

            for(int symIndex = 0; symIndex < _sa.symmetryPermutationCount(sibling->cluster2->structure); ++symIndex){
                if(_sa.symmetryTransformation(sibling->cluster2->structure, symIndex).equals(misorientation, 1e-6)){
                    _sa.clusterGraph().createClusterTransition(transition->cluster2, sibling->cluster2, misorientation, 2);
                    break;
                }
            }
        }
    }
}

void ClusterBuilder::formSuperClusters(){
    const size_t oldTransitionCount = _sa.clusterGraph().clusterTransitions().size();

    for(Cluster* cluster : _sa.clusterGraph().clusters()){
        if(!cluster || cluster->id == 0){
            continue;
        }
        cluster->rank = 0;
    }

    for(Cluster* cluster : _sa.clusterGraph().clusters()){
        if(!cluster || cluster->id == 0){
            continue;
        }
        if(cluster->structure != _context.inputCrystalType){
            processDefectCluster(cluster);
        }
    }

    const size_t newTransitionCount = _sa.clusterGraph().clusterTransitions().size();
    std::vector<ClusterTransition*> superClusterTransitions;
    superClusterTransitions.reserve(newTransitionCount - oldTransitionCount);
    for(size_t i = oldTransitionCount; i < newTransitionCount; ++i){
        superClusterTransitions.push_back(_sa.clusterGraph().clusterTransitions()[i]);
    }

    ClusterHierarchyUtils::rebuildParentHierarchy(_sa, superClusterTransitions);
}

Cluster* ClusterBuilder::startNew(int atomIndex, int structureType){
    Cluster* cluster = _sa.clusterGraph().createCluster(
        structureType,
        std::string(_sa.topologyName(structureType))
    );
    assert(cluster->id > 0);

    cluster->atomCount = 1;
    _context.atomClusters->setInt(atomIndex, cluster->id);
    return cluster;
}

Cluster* ClusterBuilder::getParentGrain(Cluster* cluster){
    return ClusterHierarchyUtils::getParentGrain(_sa, cluster);
}

}
