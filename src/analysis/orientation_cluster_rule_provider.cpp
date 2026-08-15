#include <volt/analysis/orientation_cluster_rule_provider.h>

#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/structures/cluster.h>

#include <cmath>
#include <cstdint>
#include <limits>

namespace Volt {

namespace{

constexpr double kOrientationClusterToleranceRadians = 25.0 * PI / 180.0;
constexpr double kOrientationNeighborVectorTolerance = 1e-4;
constexpr double kOrientationBondAlignmentTolerance = 0.9063077870366499;

}

OrientationClusterRuleProvider::OrientationClusterRuleProvider(
    std::shared_ptr<const std::vector<OrientationClusterAtomState>> atomStates
)
    : _atomStates(std::move(atomStates)){}

void OrientationClusterRuleProvider::initializeClusterSeed(
    const StructureAnalysis&,
    const AnalysisContext&,
    Cluster& cluster,
    int,
    int
) const{
    cluster.symmetryTransformation = 0;
}

const OrientationClusterAtomState* OrientationClusterRuleProvider::stateFor(int atomIndex) const{
    if(!_atomStates || atomIndex < 0 || atomIndex >= static_cast<int>(_atomStates->size())){
        return nullptr;
    }
    const auto& state = (*_atomStates)[static_cast<std::size_t>(atomIndex)];
    return state.valid ? &state : nullptr;
}

bool OrientationClusterRuleProvider::orthonormalizeMatrix(const Matrix3& input, Matrix3& output){
    Vector3 c0 = input.column(0);
    Vector3 c1 = input.column(1);
    Vector3 c2 = input.column(2);

    const double l0 = c0.length();
    if(l0 <= EPSILON){
        return false;
    }
    c0 /= l0;

    c1 -= c0 * c0.dot(c1);
    const double l1 = c1.length();
    if(l1 <= EPSILON){
        return false;
    }
    c1 /= l1;

    c2 -= c0 * c0.dot(c2);
    c2 -= c1 * c1.dot(c2);
    const double l2 = c2.length();
    if(l2 <= EPSILON){
        return false;
    }
    c2 /= l2;

    output = Matrix3(c0, c1, c2);
    if(output.determinant() < 0.0){
        output.column(2) = -output.column(2);
    }
    return true;
}

double OrientationClusterRuleProvider::matrixTrace(const Matrix3& matrix){
    return matrix(0, 0) + matrix(1, 1) + matrix(2, 2);
}

bool OrientationClusterRuleProvider::symmetryAllowed(std::uint64_t mask, int symmetryIndex){
    return mask == 0 || (
        symmetryIndex >= 0 &&
        symmetryIndex < 63 &&
        (mask & (std::uint64_t{1} << symmetryIndex)) != 0
    );
}

Vector3 OrientationClusterRuleProvider::transformedIdealVector(
    const StructureAnalysis& analysis,
    int structureType,
    int symmetryIndex,
    int slot
) const{
    if(slot < 0 || slot >= analysis.coordinationNumber(structureType)){
        return Vector3::Zero();
    }
    const int mappedSlot = analysis.symmetryPermutationEntry(structureType, symmetryIndex, slot);
    if(mappedSlot < 0 || mappedSlot >= analysis.coordinationNumber(structureType)){
        return Vector3::Zero();
    }
    return analysis.latticeVector(structureType, mappedSlot);
}

int OrientationClusterRuleProvider::countMatchingLocalOverlap(
    const StructureAnalysis& analysis,
    int structureType,
    int currentSymmetry,
    int neighborSlot,
    int neighborSymmetry,
    int reverseSlot
) const{
    const int coordinationNumber = analysis.coordinationNumber(structureType);
    int commonNeighborCount = 0;

    for(int currentSlot = 0; currentSlot < coordinationNumber; ++currentSlot){
        if(currentSlot == neighborSlot){
            continue;
        }

        const Vector3 expectedFromNeighbor =
            transformedIdealVector(analysis, structureType, currentSymmetry, currentSlot) -
            transformedIdealVector(analysis, structureType, currentSymmetry, neighborSlot);

        for(int neighborCandidateSlot = 0; neighborCandidateSlot < coordinationNumber; ++neighborCandidateSlot){
            if(neighborCandidateSlot == reverseSlot){
                continue;
            }
            if(!(expectedFromNeighbor - transformedIdealVector(
                analysis,
                structureType,
                neighborSymmetry,
                neighborCandidateSlot
            )).isZero(kOrientationNeighborVectorTolerance)){
                continue;
            }

            ++commonNeighborCount;
            break;
        }
    }

    return commonNeighborCount;
}

bool OrientationClusterRuleProvider::canonicalAlignmentTrace(
    const StructureAnalysis& analysis,
    int structureType,
    const OrientationClusterAtomState& currentState,
    int currentSymmetry,
    const OrientationClusterAtomState& neighborState,
    int neighborSymmetry,
    double& outTrace
) const{
    Matrix3 currentCanonical = Matrix3(
        currentState.orientation *
        analysis.symmetryTransformation(structureType, currentSymmetry).transposed()
    );
    Matrix3 neighborCanonical = Matrix3(
        neighborState.orientation *
        analysis.symmetryTransformation(structureType, neighborSymmetry).transposed()
    );

    Matrix3 relative = Matrix3(neighborCanonical.transposed() * currentCanonical);
    Matrix3 orthogonalized;
    if(!orthonormalizeMatrix(relative, orthogonalized)){
        return false;
    }

    outTrace = matrixTrace(orthogonalized);
    return true;
}

bool OrientationClusterRuleProvider::structureRequiresOrientationFallback(
    const StructureAnalysis& analysis,
    int structureType
) const{
    if(const auto it = _fallbackStructures.find(structureType); it != _fallbackStructures.end()){
        return it->second;
    }

    const int coordinationNumber = analysis.coordinationNumber(structureType);
    const int symmetryCount = std::max(1, analysis.symmetryPermutationCount(structureType));
    bool requiresFallback = false;

    for(int neighborIndex = 0; neighborIndex < coordinationNumber && !requiresFallback; ++neighborIndex){
        bool hasUsableOverlap = false;
        const Vector3 currentBond = transformedIdealVector(analysis, structureType, 0, neighborIndex);

        for(int neighborSymmetry = 0; neighborSymmetry < symmetryCount && !hasUsableOverlap; ++neighborSymmetry){
            for(int reverseSlot = 0; reverseSlot < coordinationNumber; ++reverseSlot){
                const Vector3 neighborBond = transformedIdealVector(
                    analysis,
                    structureType,
                    neighborSymmetry,
                    reverseSlot
                );
                if(!(currentBond + neighborBond).isZero(kOrientationNeighborVectorTolerance)){
                    continue;
                }
                if(countMatchingLocalOverlap(
                    analysis,
                    structureType,
                    0,
                    neighborIndex,
                    neighborSymmetry,
                    reverseSlot
                ) >= 2){
                    hasUsableOverlap = true;
                    break;
                }
            }
        }

        if(!hasUsableOverlap){
            requiresFallback = true;
        }
    }

    _fallbackStructures.emplace(structureType, requiresFallback);
    return requiresFallback;
}

ClusterRuleDecision OrientationClusterRuleProvider::tryAssignNeighbor(
    const StructureAnalysis& analysis,
    const AnalysisContext& context,
    const Cluster&,
    int currentAtomIndex,
    int neighborAtomIndex,
    int neighborIndex,
    int structureType,
    int& outNeighborSymmetry
) const{
    const auto* currentState = stateFor(currentAtomIndex);
    const auto* neighborState = stateFor(neighborAtomIndex);
    if(!currentState || !neighborState){
        return ClusterRuleDecision::Unhandled;
    }
    if(!structureRequiresOrientationFallback(analysis, structureType)){
        return ClusterRuleDecision::Unhandled;
    }
    if(context.structureTypes->getInt(neighborAtomIndex) != structureType){
        return ClusterRuleDecision::Rejected;
    }

    const int currentSymmetry = context.atomSymmetryPermutations->getInt(currentAtomIndex);
    if(currentSymmetry < 0){
        return ClusterRuleDecision::Rejected;
    }

    Vector3 spatialBond = context.simCell.wrapVector(
        context.positions->getPoint3(static_cast<std::size_t>(neighborAtomIndex)) -
        context.positions->getPoint3(static_cast<std::size_t>(currentAtomIndex))
    );
    const double spatialBondLength = spatialBond.length();
    if(spatialBondLength <= EPSILON){
        return ClusterRuleDecision::Rejected;
    }
    spatialBond /= spatialBondLength;

    const auto findBestAlignedBondSlot = [&](const OrientationClusterAtomState& state, int symmetryIndex, const Vector3& normalizedBond, double& outBestDot) {
        outBestDot = -std::numeric_limits<double>::infinity();
        int bestSlot = -1;

        const int coordinationNumber = analysis.coordinationNumber(structureType);
        for(int slot = 0; slot < coordinationNumber; ++slot){
            Vector3 predictedBond = state.orientation * transformedIdealVector(
                analysis,
                structureType,
                symmetryIndex,
                slot
            );
            const double predictedLength = predictedBond.length();
            if(predictedLength <= EPSILON){
                continue;
            }
            predictedBond /= predictedLength;

            const double dot = predictedBond.dot(normalizedBond);
            if(dot <= outBestDot){
                continue;
            }

            outBestDot = dot;
            bestSlot = slot;
        }

        return bestSlot;
    };

    double currentBondAlignment = -std::numeric_limits<double>::infinity();
    const int currentBondSlot = findBestAlignedBondSlot(
        *currentState,
        currentSymmetry,
        spatialBond,
        currentBondAlignment
    );
    if(currentBondSlot < 0 || currentBondAlignment < kOrientationBondAlignmentTolerance){
        return ClusterRuleDecision::Rejected;
    }

    const Vector3 currentBond = transformedIdealVector(analysis, structureType, currentSymmetry, currentBondSlot);
    const std::uint64_t allowedMask = static_cast<std::uint64_t>(
        context.atomAllowedSymmetryMasks->getInt64(static_cast<std::size_t>(neighborAtomIndex))
    );
    const int symmetryCount = std::max(1, analysis.symmetryPermutationCount(structureType));
    const double minTrace = 1.0 + 2.0 * std::cos(kOrientationClusterToleranceRadians);

    int bestNeighborSymmetry = -1;
    int bestOverlapCount = -1;
    double bestTrace = -std::numeric_limits<double>::infinity();
    bool bestMatchesPreferred = false;

    for(int neighborSymmetry = 0; neighborSymmetry < symmetryCount; ++neighborSymmetry){
        if(!symmetryAllowed(allowedMask, neighborSymmetry)){
            continue;
        }

        double neighborBondAlignment = -std::numeric_limits<double>::infinity();
        const int reverseSlot = findBestAlignedBondSlot(
            *neighborState,
            neighborSymmetry,
            -spatialBond,
            neighborBondAlignment
        );
        if(reverseSlot < 0 || neighborBondAlignment < kOrientationBondAlignmentTolerance){
            continue;
        }

        const Vector3 neighborBond = transformedIdealVector(
            analysis,
            structureType,
            neighborSymmetry,
            reverseSlot
        );
        if(!(currentBond + neighborBond).isZero(kOrientationNeighborVectorTolerance)){
            continue;
        }

        const int overlapCount = countMatchingLocalOverlap(
            analysis,
            structureType,
            currentSymmetry,
            currentBondSlot,
            neighborSymmetry,
            reverseSlot
        );

        double trace = -std::numeric_limits<double>::infinity();
        const bool hasOrientationTrace = canonicalAlignmentTrace(
            analysis,
            structureType,
            *currentState,
            currentSymmetry,
            *neighborState,
            neighborSymmetry,
            trace
        );

        if(!hasOrientationTrace || trace <= minTrace){
            continue;
        }
        const bool matchesPreferred =
            neighborState->preferredSymmetry >= 0 &&
            neighborState->preferredSymmetry == neighborSymmetry;

        const bool preferCandidate =
            bestNeighborSymmetry < 0 ||
            trace > bestTrace + 1e-8 ||
            (std::abs(trace - bestTrace) <= 1e-8 && matchesPreferred && !bestMatchesPreferred) ||
            (std::abs(trace - bestTrace) <= 1e-8 &&
             matchesPreferred == bestMatchesPreferred &&
             overlapCount > bestOverlapCount);

        if(!preferCandidate){
            continue;
        }

        bestNeighborSymmetry = neighborSymmetry;
        bestOverlapCount = overlapCount;
        bestTrace = trace;
        bestMatchesPreferred = matchesPreferred;
    }

    if(bestNeighborSymmetry < 0){
        return ClusterRuleDecision::Rejected;
    }

    outNeighborSymmetry = bestNeighborSymmetry;
    return ClusterRuleDecision::Accepted;
}

ClusterRuleDecision OrientationClusterRuleProvider::tryCalculateTransition(
    const StructureAnalysis& analysis,
    const AnalysisContext& context,
    int atomIndex,
    int neighborAtomIndex,
    int,
    Matrix3& outTransition
) const{
    const auto* atomState = stateFor(atomIndex);
    const auto* neighborState = stateFor(neighborAtomIndex);
    if(!atomState || !neighborState){
        return ClusterRuleDecision::Unhandled;
    }

    const int structureType = context.structureTypes->getInt(atomIndex);
    if(!structureRequiresOrientationFallback(analysis, structureType)){
        return ClusterRuleDecision::Unhandled;
    }
    if(structureType != context.structureTypes->getInt(neighborAtomIndex)){
        return ClusterRuleDecision::Unhandled;
    }

    const int atomSymmetry = context.atomSymmetryPermutations->getInt(atomIndex);
    const int neighborSymmetry = context.atomSymmetryPermutations->getInt(neighborAtomIndex);
    if(atomSymmetry < 0 || neighborSymmetry < 0){
        return ClusterRuleDecision::Rejected;
    }

    const Matrix3 atomOrientation = Matrix3(
        atomState->orientation *
        analysis.symmetryTransformation(structureType, atomSymmetry).transposed()
    );
    const Matrix3 neighborOrientation = Matrix3(
        neighborState->orientation *
        analysis.symmetryTransformation(structureType, neighborSymmetry).transposed()
    );

    Matrix3 rawTransition = Matrix3(neighborOrientation.transposed() * atomOrientation);
    if(!orthonormalizeMatrix(rawTransition, outTransition)){
        return ClusterRuleDecision::Rejected;
    }

    return ClusterRuleDecision::Accepted;
}

}
