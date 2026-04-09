#pragma once

#include <cctype>
#include <string>
#include <string_view>

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
inline constexpr CoordinationStructureType NUM_COORD_TYPES = 7;

using LatticeStructureType = int;
inline constexpr LatticeStructureType LATTICE_OTHER = 0;
inline constexpr LatticeStructureType LATTICE_SC = 1;
inline constexpr LatticeStructureType LATTICE_FCC = 2;
inline constexpr LatticeStructureType LATTICE_HCP = 3;
inline constexpr LatticeStructureType LATTICE_BCC = 4;
inline constexpr LatticeStructureType LATTICE_CUBIC_DIAMOND = 5;
inline constexpr LatticeStructureType LATTICE_HEX_DIAMOND = 6;
inline constexpr LatticeStructureType NUM_LATTICE_TYPES = 7;

enum { MAX_NEIGHBORS = 18 };

typedef unsigned int CNAPairBond;

inline const char* structureTypeName(int structure){
    switch(static_cast<StructureType>(structure)){
        case OTHER:
            return "OTHER";
        case SC:
            return "SC";
        case FCC:
            return "FCC";
        case HCP:
            return "HCP";
        case BCC:
            return "BCC";
        case CUBIC_DIAMOND:
            return "CUBIC_DIAMOND";
        case HEX_DIAMOND:
            return "HEX_DIAMOND";
        case ICO:
            return "ICO";
        case GRAPHENE:
            return "GRAPHENE";
        case CUBIC_DIAMOND_FIRST_NEIGH:
            return "CUBIC_DIAMOND_FIRST_NEIGH";
        case CUBIC_DIAMOND_SECOND_NEIGH:
            return "CUBIC_DIAMOND_SECOND_NEIGH";
        case HEX_DIAMOND_FIRST_NEIGH:
            return "HEX_DIAMOND_FIRST_NEIGH";
        case HEX_DIAMOND_SECOND_NEIGH:
            return "HEX_DIAMOND_SECOND_NEIGH";
        case L12:
            return "L12";
        case NUM_STRUCTURE_TYPES:
            break;
    }
    return "UNKNOWN";
}

inline const char* structureTypeName(StructureType structure){
    return structureTypeName(static_cast<int>(structure));
}

inline const char* latticeStructureTypeName(LatticeStructureType structure){
    switch(structure){
        case LATTICE_SC:
            return "SC";
        case LATTICE_FCC:
            return "FCC";
        case LATTICE_HCP:
            return "HCP";
        case LATTICE_BCC:
            return "BCC";
        case LATTICE_CUBIC_DIAMOND:
            return "CUBIC_DIAMOND";
        case LATTICE_HEX_DIAMOND:
            return "HEX_DIAMOND";
        case LATTICE_OTHER:
            return "OTHER";
        case NUM_LATTICE_TYPES:
            break;
    }
    return "UNKNOWN";
}

inline std::string normalizeCrystalTypeToken(std::string_view text){
    std::string normalized;
    normalized.reserve(text.size());
    for(const char character : text){
        if(character == '-' || character == ' '){
            normalized.push_back('_');
            continue;
        }
        normalized.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(character))));
    }
    return normalized;
}

inline bool parseLatticeStructureType(std::string_view text, LatticeStructureType& structure){
    const std::string normalized = normalizeCrystalTypeToken(text);
    if(normalized == "SC"){
        structure = LATTICE_SC;
        return true;
    }
    if(normalized == "FCC"){
        structure = LATTICE_FCC;
        return true;
    }
    if(normalized == "HCP"){
        structure = LATTICE_HCP;
        return true;
    }
    if(normalized == "BCC"){
        structure = LATTICE_BCC;
        return true;
    }
    if(normalized == "CUBIC_DIAMOND" || normalized == "DIAMOND"){
        structure = LATTICE_CUBIC_DIAMOND;
        return true;
    }
    if(normalized == "HEX_DIAMOND" || normalized == "HEXAGONAL_DIAMOND"){
        structure = LATTICE_HEX_DIAMOND;
        return true;
    }
    return false;
}

}
