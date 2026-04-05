#pragma once

#include <volt/analysis/structure_analysis.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace Volt::AnalysisSymmetryUtils{

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

inline std::uint64_t fullSymmetryMask(int symmetryCount){
    if(symmetryCount <= 0){
        return 0;
    }
    if(symmetryCount >= 63){
        return std::numeric_limits<std::uint64_t>::max();
    }
    return (std::uint64_t{1} << symmetryCount) - 1;
}

inline const std::vector<Matrix3>& cubicSymmetryRotations(){
    static const std::vector<Matrix3> rotations = []{
        std::vector<Matrix3> generated;
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
                    generated.push_back(rotation);
                }
            }
        }

        return generated;
    }();
    return rotations;
}

template<typename VectorContainer>
void findNonCoplanarVectors(const VectorContainer& vectors, int coordinationNumber, int indices[3], Matrix3& basis){
    basis = Matrix3::Zero();

    int found = 0;
    for(int vectorIndex = 0; vectorIndex < coordinationNumber && found < 3; ++vectorIndex){
        basis.column(found) = vectors[static_cast<std::size_t>(vectorIndex)];

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
        throw std::runtime_error("Unable to determine a non-coplanar basis for symmetry generation.");
    }
}

template<typename CanonicalVectorContainer, typename LatticeVectorContainer, typename SymmetryContainer>
void generateSymmetryPermutations(
    const CanonicalVectorContainer& canonicalVectors,
    int coordinationNumber,
    const LatticeVectorContainer& latticeVectors,
    SymmetryContainer& symmetries
){
    if(coordinationNumber == 0 || latticeVectors.empty()){
        return;
    }

    int basisIndices[3] = {-1, -1, -1};
    Matrix3 basis = Matrix3::Zero();
    findNonCoplanarVectors(canonicalVectors, coordinationNumber, basisIndices, basis);
    const Matrix3 basisInverse = basis.inverse();

    std::vector<int> permutation(latticeVectors.size());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::vector<int> lastPermutation(latticeVectors.size(), -1);

    using SymmetryType = typename SymmetryContainer::value_type;

    SymmetryType symmetry{};

    do{
        int changedFrom = static_cast<int>(
            std::mismatch(permutation.begin(), permutation.end(), lastPermutation.begin()).first - permutation.begin()
        );
        std::copy(permutation.begin(), permutation.end(), lastPermutation.begin());

        if(changedFrom <= basisIndices[2]){
            Matrix3 transformedBasis = Matrix3::Zero();
            transformedBasis.column(0) = latticeVectors[static_cast<std::size_t>(permutation[static_cast<std::size_t>(basisIndices[0])])];
            transformedBasis.column(1) = latticeVectors[static_cast<std::size_t>(permutation[static_cast<std::size_t>(basisIndices[1])])];
            transformedBasis.column(2) = latticeVectors[static_cast<std::size_t>(permutation[static_cast<std::size_t>(basisIndices[2])])];
            symmetry.transformation = transformedBasis * basisInverse;

            if(!symmetry.transformation.isOrthogonalMatrix()){
                bitmapSort(permutation.begin() + basisIndices[2] + 1, permutation.end(), static_cast<int>(permutation.size()));
                continue;
            }
            changedFrom = 0;
        }

        int sortFrom = basisIndices[2];
        int invalidFrom = changedFrom;
        for(; invalidFrom < coordinationNumber; ++invalidFrom){
            const Vector3 transformedVector = symmetry.transformation * canonicalVectors[static_cast<std::size_t>(invalidFrom)];
            if(!transformedVector.equals(latticeVectors[static_cast<std::size_t>(permutation[static_cast<std::size_t>(invalidFrom)])])){
                break;
            }
        }

        if(invalidFrom == coordinationNumber){
            symmetry.permutation.fill(-1);
            std::copy(
                permutation.begin(),
                permutation.begin() + coordinationNumber,
                symmetry.permutation.begin()
            );

            bool duplicate = false;
            for(const auto& existing : symmetries){
                if(existing.transformation.equals(symmetry.transformation)){
                    duplicate = true;
                    break;
                }
            }

            if(!duplicate){
                symmetries.push_back(symmetry);
            }
        }else{
            sortFrom = invalidFrom;
        }

        bitmapSort(permutation.begin() + sortFrom + 1, permutation.end(), static_cast<int>(permutation.size()));
    }while(std::next_permutation(permutation.begin(), permutation.end()));
}

template<typename SymmetryContainer>
void calculateSymmetryProducts(SymmetryContainer& symmetries){
    for(std::size_t s1 = 0; s1 < symmetries.size(); ++s1){
        if constexpr(requires { symmetries[s1].product; }){
            symmetries[s1].product.clear();
            symmetries[s1].product.reserve(symmetries.size());
        }

        symmetries[s1].inverseProduct.clear();
        symmetries[s1].inverseProduct.reserve(symmetries.size());

        for(std::size_t s2 = 0; s2 < symmetries.size(); ++s2){
            if constexpr(requires { symmetries[s1].product; }){
                const Matrix3 product = symmetries[s2].transformation * symmetries[s1].transformation;
                for(std::size_t candidate = 0; candidate < symmetries.size(); ++candidate){
                    if(symmetries[candidate].transformation.equals(product)){
                        symmetries[s1].product.push_back(static_cast<int>(candidate));
                        break;
                    }
                }
            }

            const Matrix3 inverseProduct =
                symmetries[s2].transformation.inverse() *
                symmetries[s1].transformation;

            for(std::size_t candidate = 0; candidate < symmetries.size(); ++candidate){
                if(symmetries[candidate].transformation.equals(inverseProduct)){
                    symmetries[s1].inverseProduct.push_back(static_cast<int>(candidate));
                    break;
                }
            }
        }
    }
}

template<typename SymmetryContainer>
int findClosestSymmetryPermutation(const SymmetryContainer& symmetries, const Matrix3& rotation){
    int bestIndex = 0;
    double bestDeviation = std::numeric_limits<double>::max();

    for(int symmetryIndex = 0; symmetryIndex < static_cast<int>(symmetries.size()); ++symmetryIndex){
        const Matrix3& symmetry = symmetries[static_cast<std::size_t>(symmetryIndex)].transformation;
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

}
