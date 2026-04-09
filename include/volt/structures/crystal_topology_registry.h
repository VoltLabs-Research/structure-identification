#pragma once

#include <volt/core/volt.h>

#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace Volt{

struct CrystalTopologySymmetry{
    Matrix3 transformation = Matrix3::Identity();
    std::vector<int> permutation;
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

}
