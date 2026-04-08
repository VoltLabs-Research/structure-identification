#include <volt/analysis/reconstructed_state_canonicalizer.h>
#include <volt/structures/crystal_topology_registry.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <vector>

namespace Volt {

bool isSymmetryAllowed(std::uint64_t mask, int symmetryIndex){
    return symmetryIndex >= 0
        && symmetryIndex < 63
        && (mask & (std::uint64_t{1} << symmetryIndex)) != 0;
}

bool symmetryAllowedWithFallback(std::uint64_t mask, int symmetryIndex){
    if(mask == 0){
        return true;
    }
    return isSymmetryAllowed(mask, symmetryIndex);
}

int selectInitialSymmetry(
    const StructureAnalysis& analysis,
    const AnalysisContext& context,
    int atomIndex,
    int structureType
){
    const int symmetryCount = analysis.symmetryPermutationCount(structureType);
    if(symmetryCount <= 0){
        return -1;
    }

    std::uint64_t allowedMask = std::uint64_t{0};
    if(context.atomAllowedSymmetryMasks){
        allowedMask = static_cast<std::uint64_t>(
            context.atomAllowedSymmetryMasks->getInt64(static_cast<std::size_t>(atomIndex))
        );
    }

    if(context.atomSymmetryPermutations){
        const int existingSymmetry = context.atomSymmetryPermutations->getInt(static_cast<std::size_t>(atomIndex));
        if(existingSymmetry >= 0 &&
           existingSymmetry < symmetryCount &&
           symmetryAllowedWithFallback(allowedMask, existingSymmetry)){
            return existingSymmetry;
        }
    }

    const int identityLikeSymmetry = analysis.findClosestSymmetryPermutation(structureType, Matrix3::Identity());
    if(identityLikeSymmetry >= 0 &&
       identityLikeSymmetry < symmetryCount &&
       symmetryAllowedWithFallback(allowedMask, identityLikeSymmetry)){
        return identityLikeSymmetry;
    }

    for(int symmetryIndex = 0; symmetryIndex < symmetryCount; ++symmetryIndex){
        if(symmetryAllowedWithFallback(allowedMask, symmetryIndex)){
            return symmetryIndex;
        }
    }
    return -1;
}

bool tryResolveNeighborSymmetry(
    const StructureAnalysis& analysis,
    const AnalysisContext& context,
    int currentAtomIndex,
    int currentSymmetry,
    int neighborAtomIndex,
    int neighborIndex,
    int structureType,
    int& outNeighborSymmetry
){
    Matrix3 currentLocalBasis;
    std::array<int, 3> matchedAtomIndices = { -1, -1, -1 };
    for(int axis = 0; axis < 3; ++axis){
        if(axis != 2){
            const int commonNeighborIndex = analysis.commonNeighborIndex(structureType, neighborIndex, axis);
            if(commonNeighborIndex < 0){
                return false;
            }
            const int matchedAtomIndex = analysis.getNeighbor(currentAtomIndex, commonNeighborIndex);
            if(matchedAtomIndex < 0){
                return false;
            }
            matchedAtomIndices[static_cast<std::size_t>(axis)] = matchedAtomIndex;

            const int commonNeighborSlot = analysis.symmetryPermutationEntry(
                structureType,
                currentSymmetry,
                commonNeighborIndex
            );
            const int bondNeighborSlot = analysis.symmetryPermutationEntry(
                structureType,
                currentSymmetry,
                neighborIndex
            );
            currentLocalBasis.column(axis) = analysis.latticeVector(structureType, commonNeighborSlot) -
                analysis.latticeVector(structureType, bondNeighborSlot);
        }else{
            matchedAtomIndices[static_cast<std::size_t>(axis)] = currentAtomIndex;
            const int bondNeighborSlot = analysis.symmetryPermutationEntry(
                structureType,
                currentSymmetry,
                neighborIndex
            );
            currentLocalBasis.column(axis) = -analysis.latticeVector(structureType, bondNeighborSlot);
        }
    }

    Matrix3 neighborCanonicalBasis;
    for(int axis = 0; axis < 3; ++axis){
        const int reverseSlot = analysis.findNeighbor(
            neighborAtomIndex,
            matchedAtomIndices[static_cast<std::size_t>(axis)]
        );
        if(reverseSlot < 0){
            return false;
        }
        neighborCanonicalBasis.column(axis) = analysis.latticeVector(structureType, reverseSlot);
    }

    Matrix3 inverseNeighborBasis;
    if(!neighborCanonicalBasis.inverse(inverseNeighborBasis)){
        return false;
    }

    const Matrix3 transition = currentLocalBasis * inverseNeighborBasis;
    std::uint64_t allowedMask = std::uint64_t{0};
    if(context.atomAllowedSymmetryMasks){
        allowedMask = static_cast<std::uint64_t>(
            context.atomAllowedSymmetryMasks->getInt64(static_cast<std::size_t>(neighborAtomIndex))
        );
    }

    const int symmetryCount = analysis.symmetryPermutationCount(structureType);
    for(int symmetryIndex = 0; symmetryIndex < symmetryCount; ++symmetryIndex){
        if(!symmetryAllowedWithFallback(allowedMask, symmetryIndex)){
            continue;
        }
        if(transition.equals(
            analysis.symmetryTransformation(structureType, symmetryIndex),
            CA_TRANSITION_MATRIX_EPSILON
        )){
            outNeighborSymmetry = symmetryIndex;
            return true;
        }
    }

    return false;
}

bool shouldVisitAtom(
    const AnalysisContext& context,
    std::size_t atomIndex,
    bool requireAssignedCluster
){
    if(!requireAssignedCluster){
        return true;
    }
    return context.atomClusters && context.atomClusters->getInt(atomIndex) != 0;
}

bool shouldTraverseNeighbor(
    const AnalysisContext& context,
    int currentClusterId,
    int neighborAtomIndex,
    bool requireAssignedCluster
){
    if(!requireAssignedCluster){
        return true;
    }
    return context.atomClusters &&
        context.atomClusters->getInt(static_cast<std::size_t>(neighborAtomIndex)) == currentClusterId;
}

double normalizedDotProduct(const Vector3& left, const Vector3& right){
    const double leftLength = left.length();
    const double rightLength = right.length();
    if(leftLength <= EPSILON || rightLength <= EPSILON){
        return -1.0;
    }
    return (left / leftLength).dot(right / rightLength);
}

bool solveMaximumWeightAssignment(
    const std::vector<std::vector<double>>& weights,
    std::vector<int>& assignment
){
    const int count = static_cast<int>(weights.size());
    if(count <= 0){
        assignment.clear();
        return true;
    }
    for(const auto& row : weights){
        if(static_cast<int>(row.size()) != count){
            return false;
        }
    }

    double maxWeight = weights.front().front();
    for(const auto& row : weights){
        for(double value : row){
            maxWeight = std::max(maxWeight, value);
        }
    }

    std::vector<double> potentialRow(static_cast<std::size_t>(count + 1), 0.0);
    std::vector<double> potentialColumn(static_cast<std::size_t>(count + 1), 0.0);
    std::vector<int> matchedColumn(static_cast<std::size_t>(count + 1), 0);
    std::vector<int> previousColumn(static_cast<std::size_t>(count + 1), 0);

    for(int row = 1; row <= count; ++row){
        matchedColumn[0] = row;
        int column0 = 0;
        std::vector<double> minReducedCost(static_cast<std::size_t>(count + 1), std::numeric_limits<double>::infinity());
        std::vector<unsigned char> used(static_cast<std::size_t>(count + 1), 0);

        do{
            used[static_cast<std::size_t>(column0)] = 1;
            const int currentRow = matchedColumn[static_cast<std::size_t>(column0)];
            double delta = std::numeric_limits<double>::infinity();
            int nextColumn = 0;

            for(int column = 1; column <= count; ++column){
                if(used[static_cast<std::size_t>(column)]){
                    continue;
                }

                const double cost = maxWeight - weights[static_cast<std::size_t>(currentRow - 1)][static_cast<std::size_t>(column - 1)];
                const double reducedCost = cost - potentialRow[static_cast<std::size_t>(currentRow)] -
                    potentialColumn[static_cast<std::size_t>(column)];
                if(reducedCost < minReducedCost[static_cast<std::size_t>(column)]){
                    minReducedCost[static_cast<std::size_t>(column)] = reducedCost;
                    previousColumn[static_cast<std::size_t>(column)] = column0;
                }
                if(minReducedCost[static_cast<std::size_t>(column)] < delta){
                    delta = minReducedCost[static_cast<std::size_t>(column)];
                    nextColumn = column;
                }
            }

            if(!std::isfinite(delta)){
                return false;
            }

            for(int column = 0; column <= count; ++column){
                if(used[static_cast<std::size_t>(column)]){
                    potentialRow[static_cast<std::size_t>(matchedColumn[static_cast<std::size_t>(column)])] += delta;
                    potentialColumn[static_cast<std::size_t>(column)] -= delta;
                }else{
                    minReducedCost[static_cast<std::size_t>(column)] -= delta;
                }
            }

            column0 = nextColumn;
        }while(matchedColumn[static_cast<std::size_t>(column0)] != 0);

        do{
            const int previous = previousColumn[static_cast<std::size_t>(column0)];
            matchedColumn[static_cast<std::size_t>(column0)] = matchedColumn[static_cast<std::size_t>(previous)];
            column0 = previous;
        }while(column0 != 0);
    }

    assignment.assign(static_cast<std::size_t>(count), -1);
    for(int column = 1; column <= count; ++column){
        const int row = matchedColumn[static_cast<std::size_t>(column)] - 1;
        if(row < 0 || row >= count){
            return false;
        }
        assignment[static_cast<std::size_t>(row)] = column - 1;
    }

    return std::find(assignment.begin(), assignment.end(), -1) == assignment.end();
}

bool remapNeighborShellByClusterGeometry(
    const StructureAnalysis& analysis,
    const AnalysisContext& context,
    std::size_t atomIndex,
    int structureType,
    int start,
    int exportableCount,
    const int* oldIndices,
    const CrystalTopologyEntry& topology,
    std::vector<int>& localSlotToExportSlot
){
    if(exportableCount <= 0 || !context.positions){
        return false;
    }

    Cluster* cluster = analysis.atomCluster(static_cast<int>(atomIndex));
    if(!cluster || cluster->structure != structureType){
        return false;
    }

    std::vector<std::vector<double>> weights(
        static_cast<std::size_t>(exportableCount),
        std::vector<double>(static_cast<std::size_t>(exportableCount), -1.0)
    );

    const Point3 center = context.positions->getPoint3(atomIndex);
    for(int localSlot = 0; localSlot < exportableCount; ++localSlot){
        const int neighborAtomIndex = oldIndices[start + localSlot];
        if(neighborAtomIndex < 0 || neighborAtomIndex >= static_cast<int>(context.atomCount())){
            return false;
        }

        const Vector3 spatialVector = context.simCell.wrapVector(
            context.positions->getPoint3(static_cast<std::size_t>(neighborAtomIndex)) - center
        );
        if(spatialVector.squaredLength() <= EPSILON){
            return false;
        }

        for(int exportSlot = 0; exportSlot < exportableCount; ++exportSlot){
            const Vector3 expectedVector = cluster->orientation *
                topology.latticeVectors[static_cast<std::size_t>(exportSlot)];
            weights[static_cast<std::size_t>(localSlot)][static_cast<std::size_t>(exportSlot)] =
                normalizedDotProduct(spatialVector, expectedVector);
        }
    }

    if(!solveMaximumWeightAssignment(weights, localSlotToExportSlot)){
        return false;
    }

    double minimumAssignedDot = std::numeric_limits<double>::infinity();
    for(int localSlot = 0; localSlot < exportableCount; ++localSlot){
        const int exportSlot = localSlotToExportSlot[static_cast<std::size_t>(localSlot)];
        if(exportSlot < 0 || exportSlot >= exportableCount){
            return false;
        }
        minimumAssignedDot = std::min(
            minimumAssignedDot,
            weights[static_cast<std::size_t>(localSlot)][static_cast<std::size_t>(exportSlot)]
        );
    }

    return minimumAssignedDot >= 0.5;
}

void canonicalizeSymmetryPermutationsImpl(
    StructureAnalysis& analysis,
    AnalysisContext& context,
    bool requireAssignedCluster
){
    if(!context.atomSymmetryPermutations || !context.structureTypes){
        return;
    }

    const std::size_t atomCount = context.atomCount();
    std::vector<int> canonicalSymmetry(atomCount, -1);

    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        if(!shouldVisitAtom(context, atomIndex, requireAssignedCluster)){
            continue;
        }
        const int structureType = context.structureTypes->getInt(atomIndex);
        if(structureType == LATTICE_OTHER){
            continue;
        }
        canonicalSymmetry[atomIndex] = -1;
    }

    std::deque<int> atomsToVisit;
    for(std::size_t seedAtomIndex = 0; seedAtomIndex < atomCount; ++seedAtomIndex){
        if(!shouldVisitAtom(context, seedAtomIndex, requireAssignedCluster)){
            continue;
        }
        const int structureType = context.structureTypes->getInt(seedAtomIndex);
        if(structureType == LATTICE_OTHER){
            continue;
        }
        if(canonicalSymmetry[seedAtomIndex] >= 0){
            continue;
        }

        const int clusterId = requireAssignedCluster && context.atomClusters
            ? context.atomClusters->getInt(seedAtomIndex)
            : 0;
        const int seedSymmetry = selectInitialSymmetry(
            analysis,
            context,
            static_cast<int>(seedAtomIndex),
            structureType
        );
        if(seedSymmetry < 0){
            continue;
        }

        canonicalSymmetry[seedAtomIndex] = seedSymmetry;
        atomsToVisit.clear();
        atomsToVisit.push_back(static_cast<int>(seedAtomIndex));

        while(!atomsToVisit.empty()){
            const int currentAtomIndex = atomsToVisit.front();
            atomsToVisit.pop_front();

            const int currentSymmetry = canonicalSymmetry[static_cast<std::size_t>(currentAtomIndex)];
            if(currentSymmetry < 0){
                continue;
            }

            const int coordinationNumber = analysis.coordinationNumber(structureType);
            for(int neighborIndex = 0; neighborIndex < coordinationNumber; ++neighborIndex){
                const int neighborAtomIndex = analysis.getNeighbor(currentAtomIndex, neighborIndex);
                if(neighborAtomIndex < 0 || neighborAtomIndex == currentAtomIndex){
                    continue;
                }
                if(!shouldTraverseNeighbor(
                    context,
                    clusterId,
                    neighborAtomIndex,
                    requireAssignedCluster
                )){
                    continue;
                }
                if(context.structureTypes->getInt(static_cast<std::size_t>(neighborAtomIndex)) != structureType){
                    continue;
                }
                if(canonicalSymmetry[static_cast<std::size_t>(neighborAtomIndex)] >= 0){
                    continue;
                }

                int resolvedSymmetry = -1;
                if(!tryResolveNeighborSymmetry(
                    analysis,
                    context,
                    currentAtomIndex,
                    currentSymmetry,
                    neighborAtomIndex,
                    neighborIndex,
                    structureType,
                    resolvedSymmetry
                )){
                    continue;
                }

                canonicalSymmetry[static_cast<std::size_t>(neighborAtomIndex)] = resolvedSymmetry;
                atomsToVisit.push_back(neighborAtomIndex);
            }
        }
    }

    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        if(!shouldVisitAtom(context, atomIndex, requireAssignedCluster)){
            continue;
        }
        const int structureType = context.structureTypes->getInt(atomIndex);
        if(structureType == LATTICE_OTHER){
            continue;
        }

        int resolvedSymmetry = canonicalSymmetry[atomIndex];
        if(resolvedSymmetry < 0){
            resolvedSymmetry = selectInitialSymmetry(
                analysis,
                context,
                static_cast<int>(atomIndex),
                structureType
            );
        }
        if(resolvedSymmetry < 0){
            continue;
        }
        context.atomSymmetryPermutations->setInt(atomIndex, resolvedSymmetry);
    }
}

void ReconstructedStateCanonicalizer::canonicalizeConnectedStructureSymmetries(
    StructureAnalysis& analysis,
    AnalysisContext& context
){
    canonicalizeSymmetryPermutationsImpl(analysis, context, false);
}

void ReconstructedStateCanonicalizer::canonicalizeNeighborShellsToExportConvention(
    StructureAnalysis& analysis,
    AnalysisContext& context
){
    if(!context.structureTypes ||
       !context.neighborCounts ||
       !context.neighborOffsets ||
       !context.neighborIndices ||
       !context.atomSymmetryPermutations){
        return;
    }

    const std::size_t atomCount = context.atomCount();
    const int* oldOffsets = context.neighborOffsets->constDataInt();
    const int* oldIndices = context.neighborIndices->constDataInt();
    const std::size_t totalNeighbors = context.neighborIndices->size();

    auto reorderedIndices = std::make_shared<ParticleProperty>(totalNeighbors, DataType::Int, 1, 0, false);
    std::copy(oldIndices, oldIndices + totalNeighbors, reorderedIndices->dataInt());

    std::vector<int> canonicalSymmetry(atomCount, -1);
    std::vector<Vector3> canonicalOverrides(
        atomCount * static_cast<std::size_t>(MAX_NEIGHBORS),
        Vector3::Zero()
    );
    const bool hasExistingOverrides = analysis.hasNeighborLatticeVectorOverrides();
    const auto& existingOverrides = analysis.neighborLatticeVectorOverrides();
    const std::size_t existingOverrideStride = analysis.neighborLatticeVectorOverrideStride();

    auto copyExistingOverrides = [&](std::size_t atomIndex, int count) {
        if(!hasExistingOverrides || existingOverrideStride == 0){
            return;
        }
        const std::size_t sourceBase = atomIndex * existingOverrideStride;
        if(sourceBase >= existingOverrides.size()){
            return;
        }
        const int copyCount = std::min(count, static_cast<int>(existingOverrideStride));
        for(int slot = 0; slot < copyCount; ++slot){
            const std::size_t sourceIndex = sourceBase + static_cast<std::size_t>(slot);
            if(sourceIndex >= existingOverrides.size()){
                break;
            }
            canonicalOverrides[
                atomIndex * static_cast<std::size_t>(MAX_NEIGHBORS) + static_cast<std::size_t>(slot)
            ] = existingOverrides[sourceIndex];
        }
    };

    auto writeFallbackVectors = [&](std::size_t atomIndex, int structureType, int count, int symmetryIndex) {
        copyExistingOverrides(atomIndex, count);

        if(structureType == LATTICE_OTHER || symmetryIndex < 0 ||
           symmetryIndex >= analysis.symmetryPermutationCount(structureType)){
            return;
        }
        const CrystalTopologyEntry* topology = crystalTopologyByStructureType(structureType);
        if(!topology || topology->coordinationNumber <= 0){
            return;
        }

        const int exportableCount = std::min(count, topology->coordinationNumber);
        for(int canonicalSlot = 0; canonicalSlot < exportableCount; ++canonicalSlot){
            const int localSlot = analysis.symmetryPermutationEntry(structureType, symmetryIndex, canonicalSlot);
            if(localSlot < 0 || localSlot >= topology->coordinationNumber){
                continue;
            }
            canonicalOverrides[
                atomIndex * static_cast<std::size_t>(MAX_NEIGHBORS) + static_cast<std::size_t>(localSlot)
            ] = topology->latticeVectors[static_cast<std::size_t>(localSlot)];
        }
    };

    int* newIndices = reorderedIndices->dataInt();
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        const int structureType = context.structureTypes->getInt(atomIndex);
        const int count = std::max(0, std::min(context.neighborCounts->getInt(atomIndex), static_cast<int>(MAX_NEIGHBORS)));
        const int start = oldOffsets[atomIndex];
        const int currentSymmetry = context.atomSymmetryPermutations->getInt(atomIndex);

        canonicalSymmetry[atomIndex] = currentSymmetry;
        if(structureType == LATTICE_OTHER || count <= 0 || start < 0){
            continue;
        }

        const CrystalTopologyEntry* topology = crystalTopologyByStructureType(structureType);
        if(!topology || topology->coordinationNumber <= 0 || topology->symmetries.empty() ||
           currentSymmetry < 0 || currentSymmetry >= analysis.symmetryPermutationCount(structureType)){
            writeFallbackVectors(atomIndex, structureType, count, currentSymmetry);
            continue;
        }

        const int exportSymmetry = topology->exportSymmetryIndex >= 0 &&
            topology->exportSymmetryIndex < static_cast<int>(topology->symmetries.size())
            ? topology->exportSymmetryIndex
            : 0;
        const int exportableCount = std::min(count, topology->coordinationNumber);

        std::vector<int> geometryAssignment;
        if(remapNeighborShellByClusterGeometry(
            analysis,
            context,
            atomIndex,
            structureType,
            start,
            exportableCount,
            oldIndices,
            *topology,
            geometryAssignment
        )){
            for(int localSlot = 0; localSlot < exportableCount; ++localSlot){
                const int exportSlot = geometryAssignment[static_cast<std::size_t>(localSlot)];
                newIndices[start + exportSlot] = oldIndices[start + localSlot];
                canonicalOverrides[
                    atomIndex * static_cast<std::size_t>(MAX_NEIGHBORS) +
                    static_cast<std::size_t>(exportSlot)
                ] = topology->latticeVectors[static_cast<std::size_t>(exportSlot)];
            }
            canonicalSymmetry[atomIndex] = exportSymmetry;
            continue;
        }

        std::vector<int> localByCanonical(static_cast<std::size_t>(topology->coordinationNumber), -1);
        for(int canonicalSlot = 0; canonicalSlot < exportableCount; ++canonicalSlot){
            const int localSlot = analysis.symmetryPermutationEntry(structureType, currentSymmetry, canonicalSlot);
            if(localSlot < 0 || localSlot >= topology->coordinationNumber){
                continue;
            }
            localByCanonical[static_cast<std::size_t>(canonicalSlot)] = localSlot;
        }

        bool valid = true;
        const std::size_t sourceOverrideBase = atomIndex * existingOverrideStride;
        for(int canonicalSlot = 0; canonicalSlot < exportableCount; ++canonicalSlot){
            const int exportLocalSlot =
                topology->symmetries[static_cast<std::size_t>(exportSymmetry)]
                    .permutation[static_cast<std::size_t>(canonicalSlot)];
            if(exportLocalSlot < 0 || exportLocalSlot >= topology->coordinationNumber){
                valid = false;
                break;
            }
            const int currentLocalSlot = localByCanonical[static_cast<std::size_t>(canonicalSlot)];
            if(currentLocalSlot < 0 || start + currentLocalSlot >= static_cast<int>(totalNeighbors)){
                valid = false;
                break;
            }

            newIndices[start + exportLocalSlot] = oldIndices[start + currentLocalSlot];
            canonicalOverrides[
                atomIndex * static_cast<std::size_t>(MAX_NEIGHBORS) + static_cast<std::size_t>(exportLocalSlot)
            ] = topology->latticeVectors[static_cast<std::size_t>(exportLocalSlot)];
        }

        if(valid){
            canonicalSymmetry[atomIndex] = exportSymmetry;
        }else{
            for(int slot = 0; slot < count; ++slot){
                newIndices[start + slot] = oldIndices[start + slot];
            }
            writeFallbackVectors(atomIndex, structureType, count, currentSymmetry);
        }
    }

    context.neighborIndices = std::move(reorderedIndices);
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        context.atomSymmetryPermutations->setInt(atomIndex, canonicalSymmetry[atomIndex]);
    }
    analysis.setNeighborLatticeVectorOverrides(
        std::move(canonicalOverrides),
        static_cast<std::size_t>(MAX_NEIGHBORS)
    );
}

void ReconstructedStateCanonicalizer::canonicalizeSymmetryPermutations(
    StructureAnalysis& analysis,
    AnalysisContext& context
){
    canonicalizeSymmetryPermutationsImpl(analysis, context, true);
}

}
