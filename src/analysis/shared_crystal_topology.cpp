#include <volt/analysis/shared_crystal_topology.h>

#include <volt/structures/lattice_vectors.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <mutex>
#include <numeric>
#include <stdexcept>

namespace Volt{

namespace{

template<typename Iterator>
void bitmapSort(Iterator begin, Iterator end, int maxValue){
    std::uint64_t bitArray = 0;
    for(Iterator it = begin; it != end; ++it){
        bitArray |= std::uint64_t{1} << (*it);
    }

    Iterator output = begin;
    for(int value = maxValue - 1; value >= 0; --value){
        if(bitArray & (std::uint64_t{1} << value)){
            *output++ = value;
        }
    }
}

int normalizeSharedStructureType(int structureType){
    switch(static_cast<StructureType>(structureType)){
        case StructureType::CUBIC_DIAMOND:
        case StructureType::CUBIC_DIAMOND_FIRST_NEIGH:
        case StructureType::CUBIC_DIAMOND_SECOND_NEIGH:
            return StructureType::CUBIC_DIAMOND;
        case StructureType::HEX_DIAMOND:
        case StructureType::HEX_DIAMOND_FIRST_NEIGH:
        case StructureType::HEX_DIAMOND_SECOND_NEIGH:
            return StructureType::HEX_DIAMOND;
        default:
            return structureType;
    }
}

std::vector<Matrix3> generateCubicSymmetryRotations(){
    std::vector<Matrix3> rotations;
    const std::array<std::array<int, 3>, 6> permutations = {{
        {0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}
    }};
    const int signs[8][3] = {
        { 1,  1,  1},
        { 1, -1, -1},
        {-1,  1, -1},
        {-1, -1,  1},
        { 1,  1, -1},
        { 1, -1,  1},
        {-1,  1,  1},
        {-1, -1, -1}
    };

    for(const auto& permutation : permutations){
        for(const auto& sign : signs){
            Matrix3 rotation;
            for(int axis = 0; axis < 3; ++axis){
                Vector3 column(0, 0, 0);
                column[permutation[axis]] = double(sign[axis]);
                rotation.column(axis) = column;
            }

            if(std::abs(rotation.determinant() - 1.0) < 1e-8){
                rotations.push_back(rotation);
            }
        }
    }

    return rotations;
}

void initializePrimitiveCell(
    SharedCrystalTopology& topology,
    const Vector3 (&primitiveCell)[3]
){
    topology.primitiveCell.column(0) = primitiveCell[0];
    topology.primitiveCell.column(1) = primitiveCell[1];
    topology.primitiveCell.column(2) = primitiveCell[2];
    topology.primitiveCellInverse = topology.primitiveCell.inverse();
}

template<typename BondPredicate>
void initializeBondedNeighbors(
    SharedCrystalTopology& topology,
    const Vector3* vectors,
    int coordinationNumber,
    BondPredicate bondPredicate
){
    topology.coordinationNumber = coordinationNumber;
    topology.commonNeighbors.fill({-1, -1});

    for(int ni1 = 0; ni1 < coordinationNumber; ++ni1){
        topology.neighborBonds.setNeighborBond(ni1, ni1, false);
        for(int ni2 = ni1 + 1; ni2 < coordinationNumber; ++ni2){
            topology.neighborBonds.setNeighborBond(ni1, ni2, bondPredicate(vectors[ni1], vectors[ni2]));
        }
    }
}

void initializeDiamondBonds(
    SharedCrystalTopology& topology,
    const Vector3* vectors,
    int coordinationNumber
){
    topology.coordinationNumber = coordinationNumber;
    topology.commonNeighbors.fill({-1, -1});

    for(int ni1 = 0; ni1 < coordinationNumber; ++ni1){
        topology.neighborBonds.setNeighborBond(ni1, ni1, false);
        const double cutoff = (ni1 < 4)
            ? (std::sqrt(3.0) * 0.25 + std::sqrt(0.5)) / 2.0
            : (1.0 + std::sqrt(0.5)) / 2.0;

        for(int ni2 = 0; ni2 < 4; ++ni2){
            if(ni1 < 4 && ni2 < 4){
                topology.neighborBonds.setNeighborBond(ni1, ni2, false);
            }
        }

        for(int ni2 = std::max(ni1 + 1, 4); ni2 < coordinationNumber; ++ni2){
            topology.neighborBonds.setNeighborBond(ni1, ni2, (vectors[ni1] - vectors[ni2]).length() < cutoff);
        }
    }
}

void findCommonNeighborsForBond(SharedCrystalTopology& topology, int neighborIndex){
    Matrix3 basis = Matrix3::Zero();
    basis.column(0) = topology.latticeVectors[static_cast<std::size_t>(neighborIndex)];

    if(topology.coordinationNumber == 6){
        for(int i1 = 0; i1 < 6; ++i1){
            if(i1 == neighborIndex || i1 == (neighborIndex ^ 1)){
                continue;
            }
            basis.column(1) = topology.latticeVectors[static_cast<std::size_t>(i1)];

            for(int i2 = i1 + 1; i2 < 6; ++i2){
                if(i2 == neighborIndex || i2 == (neighborIndex ^ 1) || i2 == (i1 ^ 1)){
                    continue;
                }
                basis.column(2) = topology.latticeVectors[static_cast<std::size_t>(i2)];
                if(std::abs(basis.determinant()) > EPSILON){
                    topology.commonNeighbors[static_cast<std::size_t>(neighborIndex)] = {i1, i2};
                    return;
                }
            }
        }
        return;
    }

    for(int i1 = 0; i1 < topology.coordinationNumber; ++i1){
        if(!topology.neighborBonds.neighborBond(neighborIndex, i1)){
            continue;
        }
        basis.column(1) = topology.latticeVectors[static_cast<std::size_t>(i1)];

        for(int i2 = i1 + 1; i2 < topology.coordinationNumber; ++i2){
            if(!topology.neighborBonds.neighborBond(neighborIndex, i2)){
                continue;
            }
            basis.column(2) = topology.latticeVectors[static_cast<std::size_t>(i2)];
            if(std::abs(basis.determinant()) > EPSILON){
                topology.commonNeighbors[static_cast<std::size_t>(neighborIndex)] = {i1, i2};
                return;
            }
        }
    }
}

void findNonCoplanarVectors(const SharedCrystalTopology& topology, int indices[3], Matrix3& basis){
    int found = 0;
    for(int vectorIndex = 0; vectorIndex < topology.coordinationNumber && found < 3; ++vectorIndex){
        basis.column(found) = topology.latticeVectors[static_cast<std::size_t>(vectorIndex)];

        if(found == 1){
            if(basis.column(0).cross(basis.column(1)).squaredLength() <= EPSILON){
                continue;
            }
        }else if(found == 2){
            if(std::abs(basis.determinant()) <= EPSILON){
                continue;
            }
        }

        indices[found++] = vectorIndex;
    }

    if(found != 3){
        throw std::runtime_error("Unable to determine a non-coplanar basis for shared crystal topology.");
    }
}

void generateGenericSymmetryPermutations(SharedCrystalTopology& topology){
    if(topology.coordinationNumber == 0 || topology.latticeVectors.empty()){
        return;
    }

    int basisIndices[3] = {-1, -1, -1};
    Matrix3 basis = Matrix3::Zero();
    findNonCoplanarVectors(topology, basisIndices, basis);
    const Matrix3 basisInverse = basis.inverse();

    std::vector<int> permutation(topology.latticeVectors.size());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::vector<int> lastPermutation(topology.latticeVectors.size(), -1);
    SharedCrystalSymmetryPermutation symmetry;

    do{
        int changedFrom = static_cast<int>(
            std::mismatch(permutation.begin(), permutation.end(), lastPermutation.begin()).first - permutation.begin()
        );
        std::copy(permutation.begin(), permutation.end(), lastPermutation.begin());

        if(changedFrom <= basisIndices[2]){
            Matrix3 transformedBasis = Matrix3::Zero();
            transformedBasis.column(0) = topology.latticeVectors[static_cast<std::size_t>(permutation[static_cast<std::size_t>(basisIndices[0])])];
            transformedBasis.column(1) = topology.latticeVectors[static_cast<std::size_t>(permutation[static_cast<std::size_t>(basisIndices[1])])];
            transformedBasis.column(2) = topology.latticeVectors[static_cast<std::size_t>(permutation[static_cast<std::size_t>(basisIndices[2])])];
            symmetry.transformation = transformedBasis * basisInverse;

            if(!symmetry.transformation.isOrthogonalMatrix()){
                bitmapSort(permutation.begin() + basisIndices[2] + 1, permutation.end(), static_cast<int>(permutation.size()));
                continue;
            }
            changedFrom = 0;
        }

        int sortFrom = basisIndices[2];
        int invalidFrom = changedFrom;
        for(; invalidFrom < topology.coordinationNumber; ++invalidFrom){
            const Vector3 transformedVector = symmetry.transformation * topology.latticeVectors[static_cast<std::size_t>(invalidFrom)];
            if(!transformedVector.equals(topology.latticeVectors[static_cast<std::size_t>(permutation[static_cast<std::size_t>(invalidFrom)])])){
                break;
            }
        }

        if(invalidFrom == topology.coordinationNumber){
            symmetry.permutation.fill(-1);
            std::copy(
                permutation.begin(),
                permutation.begin() + topology.coordinationNumber,
                symmetry.permutation.begin()
            );

            bool duplicate = false;
            for(const auto& existing : topology.symmetries){
                if(existing.transformation.equals(symmetry.transformation)){
                    duplicate = true;
                    break;
                }
            }

            if(!duplicate){
                topology.symmetries.push_back(std::move(symmetry));
            }
        }else{
            sortFrom = invalidFrom;
        }

        bitmapSort(permutation.begin() + sortFrom + 1, permutation.end(), static_cast<int>(permutation.size()));
    }while(std::next_permutation(permutation.begin(), permutation.end()));
}

void initializeSimpleCubicSymmetries(SharedCrystalTopology& topology){
    for(const Matrix3& rotation : generateCubicSymmetryRotations()){
        SharedCrystalSymmetryPermutation symmetry;
        symmetry.transformation = rotation;
        symmetry.permutation.fill(-1);

        for(int vectorIndex = 0; vectorIndex < topology.coordinationNumber; ++vectorIndex){
            const Vector3 transformedVector = rotation * topology.latticeVectors[static_cast<std::size_t>(vectorIndex)];
            for(int candidateIndex = 0; candidateIndex < topology.coordinationNumber; ++candidateIndex){
                if(transformedVector.equals(topology.latticeVectors[static_cast<std::size_t>(candidateIndex)])){
                    symmetry.permutation[static_cast<std::size_t>(vectorIndex)] = candidateIndex;
                    break;
                }
            }
        }

        topology.symmetries.push_back(std::move(symmetry));
    }
}

void calculateSymmetryProducts(SharedCrystalTopology& topology){
    for(std::size_t s1 = 0; s1 < topology.symmetries.size(); ++s1){
        topology.symmetries[s1].inverseProduct.reserve(topology.symmetries.size());
        for(std::size_t s2 = 0; s2 < topology.symmetries.size(); ++s2){
            const Matrix3 inverseProduct =
                topology.symmetries[s2].transformation.inverse() *
                topology.symmetries[s1].transformation;

            int matchIndex = 0;
            for(std::size_t candidate = 0; candidate < topology.symmetries.size(); ++candidate){
                if(topology.symmetries[candidate].transformation.equals(inverseProduct)){
                    matchIndex = static_cast<int>(candidate);
                    break;
                }
            }
            topology.symmetries[s1].inverseProduct.push_back(matchIndex);
        }
    }
}

template<typename BondPredicate>
void initializeBondedTopology(
    SharedCrystalTopology& topology,
    const Vector3* vectors,
    int coordinationNumber,
    int totalVectors,
    const Vector3 (&primitiveCell)[3],
    BondPredicate bondPredicate,
    bool useExplicitCubicSymmetries = false
){
    topology = SharedCrystalTopology{};
    topology.latticeVectors.assign(vectors, vectors + totalVectors);
    initializePrimitiveCell(topology, primitiveCell);
    initializeBondedNeighbors(topology, vectors, coordinationNumber, bondPredicate);

    for(int neighborIndex = 0; neighborIndex < coordinationNumber; ++neighborIndex){
        findCommonNeighborsForBond(topology, neighborIndex);
    }

    if(useExplicitCubicSymmetries){
        initializeSimpleCubicSymmetries(topology);
    }else{
        generateGenericSymmetryPermutations(topology);
    }
    calculateSymmetryProducts(topology);
}

void initializeDiamondTopology(
    SharedCrystalTopology& topology,
    const Vector3* vectors,
    int coordinationNumber,
    int totalVectors,
    const Vector3 (&primitiveCell)[3]
){
    topology = SharedCrystalTopology{};
    topology.latticeVectors.assign(vectors, vectors + totalVectors);
    initializePrimitiveCell(topology, primitiveCell);
    initializeDiamondBonds(topology, vectors, coordinationNumber);

    for(int neighborIndex = 0; neighborIndex < coordinationNumber; ++neighborIndex){
        findCommonNeighborsForBond(topology, neighborIndex);
    }

    generateGenericSymmetryPermutations(topology);
    calculateSymmetryProducts(topology);
}

const SharedCrystalTopology& initializeSharedTopology(int normalizedStructureType){
    static std::once_flag initFlag;
    static SharedCrystalTopology simpleCubicTopology;
    static SharedCrystalTopology faceCenteredCubicTopology;
    static SharedCrystalTopology hexagonalClosePackedTopology;
    static SharedCrystalTopology bodyCenteredCubicTopology;
    static SharedCrystalTopology cubicDiamondTopology;
    static SharedCrystalTopology hexDiamondTopology;

    std::call_once(initFlag, []() {
        initializeBondedTopology(
            simpleCubicTopology,
            SC_VECTORS,
            static_cast<int>(std::size(SC_VECTORS)),
            static_cast<int>(std::size(SC_VECTORS)),
            SC_PRIMITIVE_CELL,
            [](const Vector3&, const Vector3&) { return false; },
            true
        );
        initializeBondedTopology(
            faceCenteredCubicTopology,
            FCC_VECTORS,
            12,
            12,
            FCC_PRIMITIVE_CELL,
            [](const Vector3& v1, const Vector3& v2) {
                return (v1 - v2).length() < (std::sqrt(0.5) + 1.0) * 0.5;
            }
        );
        initializeBondedTopology(
            hexagonalClosePackedTopology,
            HCP_VECTORS,
            12,
            18,
            HCP_PRIMITIVE_CELL,
            [](const Vector3& v1, const Vector3& v2) {
                return (v1 - v2).length() < (std::sqrt(0.5) + 1.0) * 0.5;
            }
        );
        initializeBondedTopology(
            bodyCenteredCubicTopology,
            BCC_VECTORS,
            14,
            14,
            BCC_PRIMITIVE_CELL,
            [](const Vector3& v1, const Vector3& v2) {
                return (v1 - v2).length() < (1.0 + std::sqrt(2.0)) * 0.5;
            }
        );
        initializeDiamondTopology(
            cubicDiamondTopology,
            DIAMOND_CUBIC_VECTORS,
            16,
            20,
            CUBIC_DIAMOND_PRIMITIVE_CELL
        );
        initializeDiamondTopology(
            hexDiamondTopology,
            DIAMOND_HEX_VECTORS,
            16,
            32,
            HEXAGONAL_DIAMOND_PRIMITIVE_CELL
        );
    });

    switch(static_cast<StructureType>(normalizedStructureType)){
        case StructureType::SC:
            return simpleCubicTopology;
        case StructureType::FCC:
            return faceCenteredCubicTopology;
        case StructureType::HCP:
            return hexagonalClosePackedTopology;
        case StructureType::BCC:
            return bodyCenteredCubicTopology;
        case StructureType::CUBIC_DIAMOND:
            return cubicDiamondTopology;
        case StructureType::HEX_DIAMOND:
            return hexDiamondTopology;
        default:
            throw std::runtime_error("Unsupported shared crystal topology request.");
    }
}

}

const SharedCrystalTopology* sharedCrystalTopology(int structureType){
    const int normalized = normalizeSharedStructureType(structureType);
    switch(static_cast<StructureType>(normalized)){
        case StructureType::SC:
        case StructureType::FCC:
        case StructureType::HCP:
        case StructureType::BCC:
        case StructureType::CUBIC_DIAMOND:
        case StructureType::HEX_DIAMOND:
            return &initializeSharedTopology(normalized);
        default:
            return nullptr;
    }
}

int findClosestSharedCrystalSymmetryPermutation(
    const SharedCrystalTopology& topology,
    const Matrix3& rotation
){
    int bestIndex = 0;
    double bestDeviation = std::numeric_limits<double>::max();

    for(int symmetryIndex = 0; symmetryIndex < static_cast<int>(topology.symmetries.size()); ++symmetryIndex){
        const Matrix3& symmetry = topology.symmetries[static_cast<std::size_t>(symmetryIndex)].transformation;
        double deviation = 0.0;
        for(int row = 0; row < 3; ++row){
            for(int column = 0; column < 3; ++column){
                const double diff = rotation(row, column) - symmetry(row, column);
                deviation += diff * diff;
            }
        }
        if(deviation < bestDeviation){
            bestDeviation = deviation;
            bestIndex = symmetryIndex;
        }
    }

    return bestIndex;
}

bool adaptSharedCrystalTopology(
    const SharedCrystalTopology& sharedTopology,
    const std::vector<Vector3>& localVectors,
    AdaptedCrystalTopology& outTopology,
    double tolerance
){
    if(static_cast<int>(localVectors.size()) < sharedTopology.coordinationNumber){
        return false;
    }

    std::vector<int> localToCanonical(static_cast<std::size_t>(sharedTopology.coordinationNumber), -1);
    std::vector<int> canonicalToLocal(static_cast<std::size_t>(sharedTopology.coordinationNumber), -1);
    std::vector<unsigned char> usedCanonical(static_cast<std::size_t>(sharedTopology.coordinationNumber), 0);

    for(int localIndex = 0; localIndex < sharedTopology.coordinationNumber; ++localIndex){
        Vector3 localDirection = localVectors[static_cast<std::size_t>(localIndex)];
        const double localLength = localDirection.length();
        if(localLength <= EPSILON){
            return false;
        }
        localDirection /= localLength;

        int bestCanonicalIndex = -1;
        double bestError = std::numeric_limits<double>::max();
        for(int canonicalIndex = 0; canonicalIndex < sharedTopology.coordinationNumber; ++canonicalIndex){
            if(usedCanonical[static_cast<std::size_t>(canonicalIndex)]){
                continue;
            }

            Vector3 canonicalDirection = sharedTopology.latticeVectors[static_cast<std::size_t>(canonicalIndex)];
            const double canonicalLength = canonicalDirection.length();
            if(canonicalLength <= EPSILON){
                continue;
            }
            canonicalDirection /= canonicalLength;

            const double error = (localDirection - canonicalDirection).squaredLength();
            if(error < bestError){
                bestError = error;
                bestCanonicalIndex = canonicalIndex;
            }
        }

        if(bestCanonicalIndex < 0 || bestError > tolerance){
            return false;
        }

        localToCanonical[static_cast<std::size_t>(localIndex)] = bestCanonicalIndex;
        canonicalToLocal[static_cast<std::size_t>(bestCanonicalIndex)] = localIndex;
        usedCanonical[static_cast<std::size_t>(bestCanonicalIndex)] = 1;
    }

    outTopology = AdaptedCrystalTopology{};
    outTopology.coordinationNumber = sharedTopology.coordinationNumber;
    outTopology.latticeVectors.assign(
        localVectors.begin(),
        localVectors.begin() + sharedTopology.coordinationNumber
    );
    outTopology.commonNeighbors.fill({-1, -1});
    outTopology.symmetries.reserve(sharedTopology.symmetries.size());

    for(int localIndex = 0; localIndex < sharedTopology.coordinationNumber; ++localIndex){
        const int canonicalIndex = localToCanonical[static_cast<std::size_t>(localIndex)];
        for(int slot = 0; slot < 2; ++slot){
            const int canonicalNeighbor = sharedTopology.commonNeighbors[static_cast<std::size_t>(canonicalIndex)][static_cast<std::size_t>(slot)];
            outTopology.commonNeighbors[static_cast<std::size_t>(localIndex)][static_cast<std::size_t>(slot)] =
                canonicalNeighbor >= 0 ? canonicalToLocal[static_cast<std::size_t>(canonicalNeighbor)] : -1;
        }
    }

    for(const auto& sharedSymmetry : sharedTopology.symmetries){
        SharedCrystalSymmetryPermutation adaptedSymmetry;
        adaptedSymmetry.transformation = sharedSymmetry.transformation;
        adaptedSymmetry.permutation.fill(-1);

        for(int localIndex = 0; localIndex < sharedTopology.coordinationNumber; ++localIndex){
            const int canonicalIndex = localToCanonical[static_cast<std::size_t>(localIndex)];
            const int mappedCanonical = sharedSymmetry.permutation[static_cast<std::size_t>(canonicalIndex)];
            adaptedSymmetry.permutation[static_cast<std::size_t>(localIndex)] =
                mappedCanonical >= 0 ? canonicalToLocal[static_cast<std::size_t>(mappedCanonical)] : -1;
        }

        outTopology.symmetries.push_back(std::move(adaptedSymmetry));
    }

    for(std::size_t s1 = 0; s1 < outTopology.symmetries.size(); ++s1){
        outTopology.symmetries[s1].inverseProduct.clear();
        outTopology.symmetries[s1].inverseProduct.reserve(outTopology.symmetries.size());
        for(std::size_t s2 = 0; s2 < outTopology.symmetries.size(); ++s2){
            const Matrix3 inverseProduct =
                outTopology.symmetries[s2].transformation.inverse() *
                outTopology.symmetries[s1].transformation;

            int matchIndex = 0;
            for(std::size_t candidate = 0; candidate < outTopology.symmetries.size(); ++candidate){
                if(outTopology.symmetries[candidate].transformation.equals(inverseProduct)){
                    matchIndex = static_cast<int>(candidate);
                    break;
                }
            }
            outTopology.symmetries[s1].inverseProduct.push_back(matchIndex);
        }
    }

    return true;
}

}
