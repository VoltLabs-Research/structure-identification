#pragma once

#include <string_view>
#include <utility>

#include <volt/structures/crystal_topology_registry.h>

namespace Volt{

enum StructureType{
    OTHER = 0,
    SC,
    FCC,
    HCP,
    BCC,
    CUBIC_DIAMOND,
    HEX_DIAMOND,
    ICO,
    GRAPHENE,
    CUBIC_DIAMOND_FIRST_NEIGH,
    CUBIC_DIAMOND_SECOND_NEIGH,
    HEX_DIAMOND_FIRST_NEIGH,
    HEX_DIAMOND_SECOND_NEIGH,
    L12,
    NUM_STRUCTURE_TYPES
};

using CoordinationStructureType = int;
inline constexpr CoordinationStructureType COORD_OTHER = 0;
inline constexpr CoordinationStructureType COORD_SC = 1;
inline constexpr CoordinationStructureType COORD_FCC = 2;
inline constexpr CoordinationStructureType COORD_HCP = 3;
inline constexpr CoordinationStructureType COORD_BCC = 4;
inline constexpr CoordinationStructureType COORD_CUBIC_DIAMOND = 5;
inline constexpr CoordinationStructureType COORD_HEX_DIAMOND = 6;

using LatticeStructureType = int;
inline constexpr LatticeStructureType LATTICE_OTHER = 0;
inline constexpr LatticeStructureType LATTICE_SC = 1;
inline constexpr LatticeStructureType LATTICE_FCC = 2;
inline constexpr LatticeStructureType LATTICE_HCP = 3;
inline constexpr LatticeStructureType LATTICE_BCC = 4;
inline constexpr LatticeStructureType LATTICE_CUBIC_DIAMOND = 5;
inline constexpr LatticeStructureType LATTICE_HEX_DIAMOND = 6;

enum { MAX_NEIGHBORS = 16 };

typedef unsigned int CNAPairBond;

inline const char* legacyStructureTypeName(int structure){
    static constexpr std::pair<int, const char*> kLegacyNames[] = {
        {StructureType::OTHER, "OTHER"},
        {StructureType::SC, "SC"},
        {StructureType::FCC, "FCC"},
        {StructureType::HCP, "HCP"},
        {StructureType::BCC, "BCC"},
        {StructureType::CUBIC_DIAMOND, "CUBIC_DIAMOND"},
        {StructureType::HEX_DIAMOND, "HEX_DIAMOND"},
        {StructureType::ICO, "ICO"},
        {StructureType::GRAPHENE, "GRAPHENE"},
        {StructureType::CUBIC_DIAMOND_FIRST_NEIGH, "CUBIC_DIAMOND_FIRST_NEIGH"},
        {StructureType::CUBIC_DIAMOND_SECOND_NEIGH, "CUBIC_DIAMOND_SECOND_NEIGH"},
        {StructureType::HEX_DIAMOND_FIRST_NEIGH, "HEX_DIAMOND_FIRST_NEIGH"},
        {StructureType::HEX_DIAMOND_SECOND_NEIGH, "HEX_DIAMOND_SECOND_NEIGH"},
        {StructureType::L12, "L12"},
    };

    for(const auto& [identifier, name] : kLegacyNames){
        if(identifier == structure){
            return name;
        }
    }
    return "OTHER";
}

inline const char* legacyLatticeStructureTypeName(LatticeStructureType structure){
    static constexpr std::pair<LatticeStructureType, const char*> kLegacyNames[] = {
        {LATTICE_OTHER, "OTHER"},
        {LATTICE_SC, "SC"},
        {LATTICE_FCC, "FCC"},
        {LATTICE_HCP, "HCP"},
        {LATTICE_BCC, "BCC"},
        {LATTICE_CUBIC_DIAMOND, "CUBIC_DIAMOND"},
        {LATTICE_HEX_DIAMOND, "HEX_DIAMOND"},
    };

    for(const auto& [identifier, name] : kLegacyNames){
        if(identifier == structure){
            return name;
        }
    }
    return "UNKNOWN";
}

inline const char* structureTypeName(int structure){
    if(const auto* topology = crystalTopologyByStructureType(structure)){
        return topology->name.c_str();
    }
    return legacyStructureTypeName(structure);
}

inline const char* structureTypeName(StructureType structure){
    return structureTypeName(static_cast<int>(structure));
}

inline const char* latticeStructureTypeName(LatticeStructureType structure){
    if(const auto* topology = crystalTopologyByLatticeType(structure)){
        return topology->name.c_str();
    }
    return legacyLatticeStructureTypeName(structure);
}

inline bool parseLatticeStructureType(std::string_view text, LatticeStructureType& structure){
    if(const auto* topology = crystalTopologyByName(text)){
        structure = static_cast<LatticeStructureType>(topology->latticeType);
        return true;
    }
    return false;
}

}
