#pragma once

#include <string_view>

namespace Volt{

// Broad categories for crystal structures detected on an atom.
// Used to classify each atom's local arrangement into a known lattice
// type (e.g. FCC, BCC) or mark it as "OTHER" if it doesn't fit.
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

// Detailed coordination patterns for common-neighbor analysis.
// Describes how each atom's neighbors connect to each other, allowing
// distinguishing subtly different packings, like FCC vs HCP, or diamond
// vs hexagonal diamond.
// TODO: CREATE A FUNCTION STRUCTURETYPE TO COORDINATION STRUCTURE TYPE
enum CoordinationStructureType{
    // No matching coordination pattern
    COORD_OTHER = 0,
    COORD_SC,
    // 12 neighbors in 4-2-1 ring pattern
    COORD_FCC,
    // 12 neighbors in 4-2-2 ring pattern
    COORD_HCP,
    // 14 neighbors in 6-6-6 vs 4-4-4 pattern
    COORD_BCC,
    // 4 + 12 neighbors with 5-4-3 rings
    COORD_CUBIC_DIAMOND,
    // 4 + 12 neighbors with 5-4-4 rings
    COORD_HEX_DIAMOND,
    NUM_COORD_TYPES 
};

// High level lattice types for polyhedral template matching
// Defines the ideal reference lattices that the PTM algorithm
// can detect such as FCC or hexagonal diamond
enum LatticeStructureType{
    LATTICE_OTHER = 0,
    LATTICE_SC,
    LATTICE_FCC,
    LATTICE_HCP,
    LATTICE_BCC,
    // (zinc blende without basis)
    LATTICE_CUBIC_DIAMOND,
    // (wurtzite-type)
    LATTICE_HEX_DIAMOND,
    NUM_LATTICE_TYPES
};

// Maximum number of nearest neighbors supported by reconstructed analyses.
// Diamond-type structures require 16 slots to carry the 4 first-shell and 12 second-shell neighbors.
enum { MAX_NEIGHBORS = 16 };

// Bitmask type representing a bond between two common neighbors.
// Each bit in a CNAPairBond corresponds to one neighbor index;
// the union of two bits marks a neighbor-neighbor bond.
typedef unsigned int CNAPairBond;

inline const char* latticeStructureTypeName(LatticeStructureType structure){
    switch(structure){
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
        case LATTICE_SC:
            return "SC";
        case LATTICE_OTHER:
        case NUM_LATTICE_TYPES:
            break;
    }
    return "UNKNOWN";
}

inline bool parseLatticeStructureType(std::string_view text, LatticeStructureType& structure){
    if(text == "FCC"){
        structure = LATTICE_FCC;
        return true;
    }
    if(text == "BCC"){
        structure = LATTICE_BCC;
        return true;
    }
    if(text == "HCP"){
        structure = LATTICE_HCP;
        return true;
    }
    if(text == "SC"){
        structure = LATTICE_SC;
        return true;
    }
    if(text == "CUBIC_DIAMOND"){
        structure = LATTICE_CUBIC_DIAMOND;
        return true;
    }
    if(text == "HEX_DIAMOND"){
        structure = LATTICE_HEX_DIAMOND;
        return true;
    }
    return false;
}

}
