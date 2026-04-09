#include <volt/structures/crystal_topology_registry.h>

#include <volt/analysis/crystal_symmetry_utils.h>

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <system_error>

#if defined(__linux__)
#include <unistd.h>
#endif

namespace Volt{

namespace{

std::mutex& topologyRegistryMutex(){
    static std::mutex mutex;
    return mutex;
}

std::vector<std::filesystem::path>& configuredTopologyRoots(){
    static std::vector<std::filesystem::path> roots;
    return roots;
}

std::string normalizeKey(std::string value){
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char character){
        if(character == '-' || character == ' '){
            return '_';
        }
        return static_cast<char>(std::tolower(character));
    });
    return value;
}

Vector3 parseVector3(const YAML::Node& value){
    if(!value || !value.IsSequence() || value.size() != 3){
        throw std::runtime_error("Expected a 3-component vector.");
    }
    return Vector3(
        value[0].as<double>(),
        value[1].as<double>(),
        value[2].as<double>()
    );
}

std::vector<Vector3> parseVectorList(const YAML::Node& value){
    std::vector<Vector3> result;
    if(!value || !value.IsSequence()){
        return result;
    }
    result.reserve(value.size());
    for(const auto& item : value){
        result.push_back(parseVector3(item));
    }
    return result;
}

std::vector<Vector3> activeNeighborShell(const CrystalTopologyEntry& entry){
    return std::vector<Vector3>(
        entry.neighborVectors.begin(),
        entry.neighborVectors.begin() + entry.coordinationNumber
    );
}

std::vector<CrystalTopologySymmetry> generateSymmetryPermutations(const CrystalTopologyEntry& entry){
    std::vector<CrystalTopologySymmetry> symmetries;
    const std::vector<Vector3> canonicalNeighborVectors = activeNeighborShell(entry);
    AnalysisSymmetryUtils::generateSymmetryPermutations(
        canonicalNeighborVectors,
        entry.coordinationNumber,
        canonicalNeighborVectors,
        symmetries
    );
    AnalysisSymmetryUtils::retainProperRotations(symmetries);
    return symmetries;
}

std::vector<std::filesystem::path> normalizeSearchRoots(std::vector<std::filesystem::path> roots){
    std::vector<std::filesystem::path> normalized;
    normalized.reserve(roots.size());
    for(auto& root : roots){
        if(root.empty()){
            continue;
        }
        const auto canonical = root.lexically_normal();
        if(std::find(normalized.begin(), normalized.end(), canonical) == normalized.end()){
            normalized.push_back(canonical);
        }
    }
    return normalized;
}

std::vector<std::filesystem::path> defaultTopologySearchRoots(){
    if(const char* envRoot = std::getenv("VOLT_LATTICE_DIR")){
        if(*envRoot != '\0'){
            return normalizeSearchRoots({std::filesystem::path(envRoot)});
        }
    }

    std::vector<std::filesystem::path> roots;

#ifdef VOLT_OPENDXA_LATTICE_SOURCE_DIR
    roots.emplace_back(std::filesystem::path(VOLT_OPENDXA_LATTICE_SOURCE_DIR));
#endif

#if defined(__linux__)
    std::array<char, 4096> buffer{};
    const ssize_t length = ::readlink("/proc/self/exe", buffer.data(), buffer.size() - 1);
    if(length > 0){
        buffer[static_cast<std::size_t>(length)] = '\0';
        std::filesystem::path executablePath(buffer.data());
        std::filesystem::path current = executablePath.parent_path();
        for(int depth = 0; depth < 6 && !current.empty(); ++depth){
            roots.push_back(current / "share/volt/lattices");
            current = current.parent_path();
        }
    }
#endif

    return normalizeSearchRoots(std::move(roots));
}

std::vector<std::filesystem::path> topologySearchRootsUnsafe(){
    const auto configured = configuredTopologyRoots();
    if(!configured.empty()){
        return configured;
    }
    return defaultTopologySearchRoots();
}

std::vector<std::filesystem::path> discoverTopologyFiles(const std::vector<std::filesystem::path>& searchRoots){
    std::vector<std::filesystem::path> topologyFiles;
    std::error_code error;

    for(const auto& root : searchRoots){
        if(!std::filesystem::exists(root, error) || !std::filesystem::is_directory(root, error)){
            continue;
        }

        for(const auto& entry : std::filesystem::recursive_directory_iterator(root, error)){
            if(error){
                break;
            }
            if(!entry.is_regular_file()){
                continue;
            }
            const auto extension = normalizeKey(entry.path().extension().string());
            if(extension != ".yml" && extension != ".yaml"){
                continue;
            }
            topologyFiles.push_back(entry.path());
        }
    }

    std::sort(topologyFiles.begin(), topologyFiles.end());
    topologyFiles.erase(std::unique(topologyFiles.begin(), topologyFiles.end()), topologyFiles.end());
    return topologyFiles;
}

CrystalTopologyEntry parseTopologyEntry(const std::filesystem::path& filePath){
    YAML::Node document;
    try{
        document = YAML::LoadFile(filePath.string());
    }catch(const std::exception& error){
        throw std::runtime_error(
            "Unable to parse topology YAML '" + filePath.string() + "': " + error.what()
        );
    }

    if(!document || !document.IsMap()){
        throw std::runtime_error("Topology YAML must contain a mapping root.");
    }

    CrystalTopologyEntry entry;
    entry.name = document["name"].as<std::string>();
    entry.coordinationNumber = document["coordination_number"].as<int>();
    entry.neighborVectors = parseVectorList(document["neighbor_vectors"]);

    if(entry.name.empty()){
        throw std::runtime_error("Topology name cannot be empty.");
    }
    if(entry.coordinationNumber <= 0){
        throw std::runtime_error("coordination_number must be positive.");
    }
    if(static_cast<int>(entry.neighborVectors.size()) < entry.coordinationNumber){
        throw std::runtime_error("neighbor_vectors count is smaller than coordination_number.");
    }
    entry.symmetries = generateSymmetryPermutations(entry);
    if(entry.symmetries.empty()){
        throw std::runtime_error("Failed to derive symmetry permutations from neighbor_vectors.");
    }

    return entry;
}

void registerTopologyName(std::unordered_map<std::string, std::size_t>& nameIndex, const std::string& name, std::size_t index){
    if(name.empty()){
        return;
    }
    const auto normalized = normalizeKey(name);
    const auto [it, inserted] = nameIndex.emplace(normalized, index);
    if(!inserted && it->second != index){
        throw std::runtime_error("Duplicate topology alias: " + name);
    }
}

bool almostEqual(double lhs, double rhs){
    return std::abs(lhs - rhs) <= 1e-12;
}

bool equivalentVector3(const Vector3& lhs, const Vector3& rhs){
    return almostEqual(lhs.x(), rhs.x())
        && almostEqual(lhs.y(), rhs.y())
        && almostEqual(lhs.z(), rhs.z());
}

bool equivalentMatrix3(const Matrix3& lhs, const Matrix3& rhs){
    for(int row = 0; row < 3; ++row){
        for(int column = 0; column < 3; ++column){
            if(!almostEqual(lhs(row, column), rhs(row, column))){
                return false;
            }
        }
    }
    return true;
}

bool equivalentTopologyEntries(const CrystalTopologyEntry& lhs, const CrystalTopologyEntry& rhs){
    if(lhs.name != rhs.name ||
       lhs.coordinationNumber != rhs.coordinationNumber ||
       lhs.neighborVectors.size() != rhs.neighborVectors.size() ||
       lhs.symmetries.size() != rhs.symmetries.size()){
        return false;
    }

    for(std::size_t index = 0; index < lhs.neighborVectors.size(); ++index){
        if(!equivalentVector3(lhs.neighborVectors[index], rhs.neighborVectors[index])){
            return false;
        }
    }
    for(std::size_t index = 0; index < lhs.symmetries.size(); ++index){
        if(!equivalentMatrix3(lhs.symmetries[index].transformation, rhs.symmetries[index].transformation) ||
           lhs.symmetries[index].permutation != rhs.symmetries[index].permutation){
            return false;
        }
    }

    return true;
}

}

CrystalTopologyRegistry::CrystalTopologyRegistry(std::vector<std::filesystem::path> searchRoots){
    const auto topologyFiles = discoverTopologyFiles(searchRoots);
    if(topologyFiles.empty()){
        throw std::runtime_error("No topology YAML files were found.");
    }

    _entries.reserve(topologyFiles.size());
    std::unordered_map<std::string, std::size_t> existingByName;
    for(const auto& filePath : topologyFiles){
        CrystalTopologyEntry entry = parseTopologyEntry(filePath);
        const std::string normalizedName = normalizeKey(entry.name);
        const auto existing = existingByName.find(normalizedName);
        if(existing == existingByName.end()){
            existingByName.emplace(normalizedName, _entries.size());
            _entries.push_back(std::move(entry));
            continue;
        }

        if(!equivalentTopologyEntries(_entries[existing->second], entry)){
            throw std::runtime_error(
                "Conflicting topology definitions found for lattice '" + entry.name + "'."
            );
        }
    }

    for(std::size_t index = 0; index < _entries.size(); ++index){
        const auto& entry = _entries[index];
        registerTopologyName(_nameIndex, entry.name, index);
    }
}

const CrystalTopologyRegistry& CrystalTopologyRegistry::instance(){
    static std::unique_ptr<CrystalTopologyRegistry> registry;
    static std::vector<std::filesystem::path> loadedRoots;

    std::lock_guard<std::mutex> lock(topologyRegistryMutex());
    const auto searchRoots = topologySearchRootsUnsafe();
    if(!registry || loadedRoots != searchRoots){
        registry.reset(new CrystalTopologyRegistry(searchRoots));
        loadedRoots = searchRoots;
    }
    return *registry;
}

const CrystalTopologyEntry* CrystalTopologyRegistry::findByName(std::string_view name) const{
    const auto it = _nameIndex.find(normalizeKey(std::string(name)));
    return it == _nameIndex.end() ? nullptr : &_entries[it->second];
}

const CrystalTopologyRegistry& crystalTopologyRegistry(){
    return CrystalTopologyRegistry::instance();
}

void setCrystalTopologySearchRoots(std::vector<std::filesystem::path> roots){
    std::lock_guard<std::mutex> lock(topologyRegistryMutex());
    configuredTopologyRoots() = normalizeSearchRoots(std::move(roots));
}

void setCrystalTopologySearchRoot(std::filesystem::path root){
    setCrystalTopologySearchRoots({std::move(root)});
}

const CrystalTopologyEntry* crystalTopologyByName(std::string_view name){
    return crystalTopologyRegistry().findByName(name);
}

}
