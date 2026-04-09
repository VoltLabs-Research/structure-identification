#include <volt/analysis/crystal_topology_library.h>
#include <volt/analysis/crystal_symmetry_utils.h>
#include <volt/structures/crystal_topology_registry.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <unordered_map>

namespace Volt{

namespace SharedCrystalTopologyDetail{

const CrystalTopologyEntry& requireTopologyEntry(std::string_view topologyName){
    const auto* entry = crystalTopologyByName(topologyName);
    if(!entry){
        throw std::runtime_error("Unsupported shared crystal topology request.");
    }
    return *entry;
}

void initializeCommonNeighbors(SharedCrystalTopology& topology){
    topology.commonNeighbors.fill({-1, -1});

    for(int neighborIndex = 0; neighborIndex < topology.coordinationNumber; ++neighborIndex){
        Matrix3 basis = Matrix3::Zero();
        basis.column(0) = topology.neighborVectors[static_cast<std::size_t>(neighborIndex)];

        double minBondDistanceSq = std::numeric_limits<double>::max();
        for(int otherIndex = 0; otherIndex < topology.coordinationNumber; ++otherIndex){
            if(otherIndex == neighborIndex){
                continue;
            }
            const double distSq = (
                topology.neighborVectors[static_cast<std::size_t>(neighborIndex)] -
                topology.neighborVectors[static_cast<std::size_t>(otherIndex)]
            ).squaredLength();
            if(distSq > EPSILON && distSq < minBondDistanceSq){
                minBondDistanceSq = distSq;
            }
        }

        if(!std::isfinite(minBondDistanceSq)){
            continue;
        }

        for(int i1 = 0; i1 < topology.coordinationNumber; ++i1){
            if(i1 == neighborIndex){
                continue;
            }
            const double d1 = (
                topology.neighborVectors[static_cast<std::size_t>(neighborIndex)] -
                topology.neighborVectors[static_cast<std::size_t>(i1)]
            ).squaredLength();
            if(std::abs(d1 - minBondDistanceSq) > 1e-6){
                continue;
            }
            basis.column(1) = topology.neighborVectors[static_cast<std::size_t>(i1)];

            for(int i2 = i1 + 1; i2 < topology.coordinationNumber; ++i2){
                if(i2 == neighborIndex){
                    continue;
                }
                const double d2 = (
                    topology.neighborVectors[static_cast<std::size_t>(neighborIndex)] -
                    topology.neighborVectors[static_cast<std::size_t>(i2)]
                ).squaredLength();
                if(std::abs(d2 - minBondDistanceSq) > 1e-6){
                    continue;
                }
                basis.column(2) = topology.neighborVectors[static_cast<std::size_t>(i2)];
                if(std::abs(basis.determinant()) > EPSILON){
                    topology.commonNeighbors[static_cast<std::size_t>(neighborIndex)] = {i1, i2};
                    goto next_neighbor;
                }
            }
        }

    next_neighbor:
        continue;
    }
}

void initializeExplicitSymmetries(SharedCrystalTopology& topology, const CrystalTopologyEntry& entry){
    topology.symmetries.clear();
    topology.symmetries.reserve(entry.symmetries.size());

    for(const CrystalTopologySymmetry& source : entry.symmetries){
        SharedCrystalSymmetryPermutation symmetry;
        symmetry.transformation = source.transformation;
        symmetry.permutation.fill(-1);

        for(int index = 0; index < topology.coordinationNumber; ++index){
            symmetry.permutation[static_cast<std::size_t>(index)] =
                source.permutation[static_cast<std::size_t>(index)];
        }

        topology.symmetries.push_back(std::move(symmetry));
    }
}

void initializeTopologyFromEntry(SharedCrystalTopology& topology, const CrystalTopologyEntry& entry){
    topology = SharedCrystalTopology{};
    topology.coordinationNumber = entry.coordinationNumber;
    topology.neighborVectors = entry.neighborVectors;
    initializeCommonNeighbors(topology);
    initializeExplicitSymmetries(topology, entry);
    AnalysisSymmetryUtils::calculateSymmetryProducts(topology.symmetries);
}

const SharedCrystalTopology& initializeSharedTopology(std::string_view topologyName){
    static std::once_flag initFlag;
    static std::unordered_map<std::string, SharedCrystalTopology> topologies;

    std::call_once(initFlag, []() {
        for(const CrystalTopologyEntry& entry : crystalTopologyRegistry().entries()){
            initializeTopologyFromEntry(topologies[entry.name], entry);
        }
    });

    const CrystalTopologyEntry& entry = requireTopologyEntry(topologyName);
    const auto it = topologies.find(entry.name);
    if(it == topologies.end()){
        throw std::runtime_error("Unsupported shared crystal topology request.");
    }
    return it->second;
}

}

using namespace SharedCrystalTopologyDetail;

const SharedCrystalTopology* sharedCrystalTopology(std::string_view topologyName){
    return crystalTopologyByName(topologyName)
        ? &initializeSharedTopology(topologyName)
        : nullptr;
}

int findClosestSharedCrystalSymmetryPermutation(
    const SharedCrystalTopology& topology,
    const Matrix3& rotation
){
    return AnalysisSymmetryUtils::findClosestSymmetryPermutation(topology.symmetries, rotation);
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

            Vector3 canonicalDirection = sharedTopology.neighborVectors[static_cast<std::size_t>(canonicalIndex)];
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
    outTopology.neighborVectors.assign(
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

    AnalysisSymmetryUtils::calculateSymmetryProducts(outTopology.symmetries);
    return true;
}

}
