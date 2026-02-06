#pragma once

#include <vector>
#include <array>
#include <volt/math/lin_alg.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/structures/coordination_structure.h>

namespace Volt {

struct SymmetryPermutation{
    Matrix3 transformation;
    std::array<int, MAX_NEIGHBORS> permutation;
    std::vector<int> product;
    std::vector<int> inverseProduct;
};

struct LatticeStructure{
    const CoordinationStructure* coordStructure;
    std::vector<Vector3> latticeVectors;
    Matrix3 primitiveCell;
    Matrix3 primitiveCellInverse;
    int maxNeighbors;
    std::vector<SymmetryPermutation> permutations;
};

}
