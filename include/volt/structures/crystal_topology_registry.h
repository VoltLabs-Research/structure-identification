#pragma once

#include <volt/core/volt.h>
#include <volt/structures/crystal_structure_types.h>

#include <array>
#include <cstddef>
#include <filesystem>
#include <mutex>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace YAML{
class Node;
}

namespace Volt{

struct CrystalTopologySymmetry{
    Matrix3 transformation = Matrix3::Identity();
    std::array<int, MAX_NEIGHBORS> permutation{};
};

struct CrystalTopologyEntry{
    std::string name;
    int coordinationNumber = 0;
    std::vector<Vector3> neighborVectors;
    std::vector<CrystalTopologySymmetry> symmetries;
};

class CrystalTopologyRegistry{
public:
    static const CrystalTopologyRegistry& instance();

    const CrystalTopologyEntry* findByName(std::string_view name) const;

    const std::vector<CrystalTopologyEntry>& entries() const{
        return _entries;
    }

private:
    explicit CrystalTopologyRegistry(std::vector<std::filesystem::path> searchRoots);

    std::vector<CrystalTopologyEntry> _entries;
    std::unordered_map<std::string, std::size_t> _nameIndex;
};

void setCrystalTopologySearchRoots(std::vector<std::filesystem::path> roots);
void setCrystalTopologySearchRoot(std::filesystem::path root);
const CrystalTopologyRegistry& crystalTopologyRegistry();
const CrystalTopologyEntry* crystalTopologyByName(std::string_view name);

// What the registry is built from: the search roots, the topology files and their entries.
namespace CrystalTopologyRegistryDetail{

std::mutex& topologyRegistryMutex();
std::vector<std::filesystem::path>& configuredTopologyRoots();
std::string normalizeKey(std::string value);
Vector3 parseVector3(const YAML::Node& value);
std::vector<Vector3> parseVectorList(const YAML::Node& value);
std::vector<Vector3> activeNeighborShell(const CrystalTopologyEntry& entry);
std::vector<CrystalTopologySymmetry> generateSymmetryPermutations(const CrystalTopologyEntry& entry);
std::vector<std::filesystem::path> normalizeSearchRoots(std::vector<std::filesystem::path> roots);
std::vector<std::filesystem::path> defaultTopologySearchRoots();
std::vector<std::filesystem::path> topologySearchRootsUnsafe();
std::vector<std::filesystem::path> discoverTopologyFiles(const std::vector<std::filesystem::path>& searchRoots);
CrystalTopologyEntry parseTopologyEntry(const std::filesystem::path& filePath);
void registerTopologyName(std::unordered_map<std::string, std::size_t>& nameIndex, const std::string& name, std::size_t index);
bool almostEqual(double lhs, double rhs);
bool equivalentVector3(const Vector3& lhs, const Vector3& rhs);
bool equivalentMatrix3(const Matrix3& lhs, const Matrix3& rhs);
bool equivalentTopologyEntries(const CrystalTopologyEntry& lhs, const CrystalTopologyEntry& rhs);

}

}
