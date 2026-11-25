#pragma once

#include <volt/core/volt.h>
#include <volt/structures/crystal_structure_types.h>

#include <array>
#include <string_view>
#include <vector>

namespace Volt{

struct CrystalTopologyEntry;

struct SharedCrystalSymmetryPermutation{
    Matrix3 transformation = Matrix3::Identity();
    std::array<int, MAX_NEIGHBORS> permutation{};
    std::vector<int> inverseProduct;
};

struct SharedCrystalTopology{
    int coordinationNumber = 0;
    std::vector<Vector3> neighborVectors;
    std::array<std::array<int, 2>, MAX_NEIGHBORS> commonNeighbors{};
    std::vector<SharedCrystalSymmetryPermutation> symmetries;
};

struct AdaptedCrystalTopology{
    int coordinationNumber = 0;
    std::vector<Vector3> neighborVectors;
    std::array<std::array<int, 2>, MAX_NEIGHBORS> commonNeighbors{};
    std::vector<SharedCrystalSymmetryPermutation> symmetries;
};

const SharedCrystalTopology* sharedCrystalTopology(std::string_view topologyName);
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

// How a shared topology is built from its registry entry.
namespace SharedCrystalTopologyDetail{

const CrystalTopologyEntry& requireTopologyEntry(std::string_view topologyName);
void initializeTopologyFromEntry(SharedCrystalTopology& topology, const CrystalTopologyEntry& entry);
void initializeCommonNeighbors(SharedCrystalTopology& topology);
void initializeExplicitSymmetries(SharedCrystalTopology& topology, const CrystalTopologyEntry& entry);
const SharedCrystalTopology& initializeSharedTopology(std::string_view topologyName);

}

}
