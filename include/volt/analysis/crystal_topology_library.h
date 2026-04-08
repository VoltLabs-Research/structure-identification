#pragma once

#include <volt/core/volt.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/structures/neighbor_bond_array.h>

#include <array>
#include <vector>

namespace Volt{

struct SharedCrystalSymmetryPermutation{
    Matrix3 transformation = Matrix3::Identity();
    std::array<int, MAX_NEIGHBORS> permutation{};
    std::vector<int> inverseProduct;
};

struct SharedCrystalTopology{
    int coordinationNumber = 0;
    std::vector<Vector3> latticeVectors;
    NeighborBondArray neighborBonds;
    std::array<std::array<int, 2>, MAX_NEIGHBORS> commonNeighbors{};
    Matrix3 primitiveCell = Matrix3::Zero();
    Matrix3 primitiveCellInverse = Matrix3::Zero();
    std::vector<SharedCrystalSymmetryPermutation> symmetries;
};

struct AdaptedCrystalTopology{
    int coordinationNumber = 0;
    std::vector<Vector3> latticeVectors;
    std::array<std::array<int, 2>, MAX_NEIGHBORS> commonNeighbors{};
    std::vector<SharedCrystalSymmetryPermutation> symmetries;
    std::vector<int> localToCanonical;
    std::vector<int> canonicalToLocal;
};

const SharedCrystalTopology* sharedCrystalTopology(int structureType);
int findClosestSharedCrystalSymmetryPermutation(
    const SharedCrystalTopology& topology,
    const Matrix3& rotation
);
bool adaptSharedCrystalTopology(
    const SharedCrystalTopology& sharedTopology,
    const std::vector<Vector3>& localVectors,
    AdaptedCrystalTopology& outTopology,
    double tolerance = 1e-4
);

}
