#include <volt/analysis/crystal_topology_library.h>
#include <volt/analysis/crystal_symmetry_utils.h>

#include <volt/structures/crystal_topology_registry.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

namespace Volt{

namespace SharedCrystalTopologyDetail{

const CrystalTopologyEntry& requireTopologyEntry(int structureType){
    const auto* entry = crystalTopologyByStructureType(structureType);
    if(!entry){
        throw std::runtime_error("Unsupported shared crystal topology request.");
    }
    return *entry;
}

void initializePrimitiveCell(SharedCrystalTopology& topology, const CrystalTopologyEntry& entry){
    topology.primitiveCell = entry.primitiveCell;
    topology.primitiveCellInverse = entry.primitiveCellInverse;
}

void initializeBondedNeighbors(SharedCrystalTopology& topology, const CrystalTopologyEntry& entry){
    topology.coordinationNumber = entry.coordinationNumber;
    for(int ni1 = 0; ni1 < entry.coordinationNumber; ++ni1){
        topology.neighborBonds.neighborArray[static_cast<std::size_t>(ni1)] =
            entry.neighborBondRows[static_cast<std::size_t>(ni1)];
        topology.commonNeighbors[static_cast<std::size_t>(ni1)] =
            entry.commonNeighbors[static_cast<std::size_t>(ni1)];
    }
}

void initializeExplicitSymmetries(SharedCrystalTopology& topology, const CrystalTopologyEntry& entry){
    topology.symmetries.clear();
    topology.symmetries.reserve(entry.symmetries.size());
    for(const auto& symmetry : entry.symmetries){
        SharedCrystalSymmetryPermutation sharedSymmetry;
        sharedSymmetry.transformation = symmetry.transformation;
        sharedSymmetry.permutation.fill(-1);
        for(std::size_t slot = 0; slot < symmetry.permutation.size() && slot < sharedSymmetry.permutation.size(); ++slot){
            sharedSymmetry.permutation[slot] = symmetry.permutation[slot];
        }
        topology.symmetries.push_back(std::move(sharedSymmetry));
    }
}

void calculateSymmetryProducts(SharedCrystalTopology& topology){
    AnalysisSymmetryUtils::calculateSymmetryProducts(topology.symmetries);
}

void initializeTopologyFromEntry(SharedCrystalTopology& topology, const CrystalTopologyEntry& entry){
    topology = SharedCrystalTopology{};
    topology.latticeVectors = entry.latticeVectors;
    initializePrimitiveCell(topology, entry);
    initializeBondedNeighbors(topology, entry);
    initializeExplicitSymmetries(topology, entry);
    calculateSymmetryProducts(topology);
}

const SharedCrystalTopology& initializeSharedTopology(int structureType){
    static std::once_flag initFlag;
    static std::unordered_map<int, SharedCrystalTopology> topologies;

    std::call_once(initFlag, []() {
        for(const CrystalTopologyEntry& entry : crystalTopologyRegistry().entries()){
            if(entry.structureType <= 0){
                continue;
            }
            initializeTopologyFromEntry(topologies[entry.structureType], entry);
        }
    });

    const CrystalTopologyEntry& entry = requireTopologyEntry(structureType);
    const auto it = topologies.find(entry.structureType);
    if(it == topologies.end()){
        throw std::runtime_error("Unsupported shared crystal topology request.");
    }
    return it->second;
}

}

using namespace SharedCrystalTopologyDetail;

const SharedCrystalTopology* sharedCrystalTopology(int structureType){
    return crystalTopologyByStructureType(structureType)
        ? &initializeSharedTopology(structureType)
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

    Matrix3 canonicalBasis = Matrix3::Zero();
    Matrix3 localBasis = Matrix3::Zero();
    int basisCount = 0;
    for(int localIndex = 0; localIndex < sharedTopology.coordinationNumber && basisCount < 3; ++localIndex){
        const int canonicalIndex = localToCanonical[static_cast<std::size_t>(localIndex)];
        if(canonicalIndex < 0){
            return false;
        }

        canonicalBasis.column(basisCount) = sharedTopology.latticeVectors[static_cast<std::size_t>(canonicalIndex)];
        localBasis.column(basisCount) = localVectors[static_cast<std::size_t>(localIndex)];

        if(basisCount == 1){
            if(canonicalBasis.column(0).cross(canonicalBasis.column(1)).squaredLength() <= EPSILON){
                continue;
            }
            if(localBasis.column(0).cross(localBasis.column(1)).squaredLength() <= EPSILON){
                continue;
            }
        }else if(basisCount == 2){
            if(std::abs(canonicalBasis.determinant()) <= EPSILON){
                continue;
            }
            if(std::abs(localBasis.determinant()) <= EPSILON){
                continue;
            }
        }

        ++basisCount;
    }

    if(basisCount != 3){
        return false;
    }

    Matrix3 canonicalBasisInverse;
    if(!canonicalBasis.inverse(canonicalBasisInverse)){
        return false;
    }

    const Matrix3 frameCanonicalToLocal = localBasis * canonicalBasisInverse;
    Matrix3 frameLocalToCanonical;
    if(!frameCanonicalToLocal.inverse(frameLocalToCanonical)){
        return false;
    }

    outTopology = AdaptedCrystalTopology{};
    outTopology.coordinationNumber = sharedTopology.coordinationNumber;
    outTopology.localToCanonical = localToCanonical;
    outTopology.canonicalToLocal = canonicalToLocal;
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
        adaptedSymmetry.transformation =
            frameCanonicalToLocal * sharedSymmetry.transformation * frameLocalToCanonical;
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
