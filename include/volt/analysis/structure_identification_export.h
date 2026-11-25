#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include <volt/utilities/parquet_atom_writer.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/core/lammps_parser.h>

namespace Volt::StructureIdentificationExport {

using json = nlohmann::json;

// The name of a structure type in the atoms table when no resolver names it.
std::string defaultStructureName(int structureType);
// The name of a structure type in the tables of an analysis built on the structures (grains, strain): OTHER for one
// out of range and for one no structure identification of theirs gives (L12).
const char* structureTypeNameForExport(int structureType);

// The atoms of each structure type, the types that have any sorted by export name, and how many atoms are in a
// cluster (or a grain): what the per-structure atoms table lists.
struct AtomsByStructure{
    std::vector<std::vector<std::size_t>> atomsOfType;
    std::vector<int> order;
    int clusteredAtoms = 0;
    std::size_t atoms = 0;
};

// A type out of range counts as OTHER.
AtomsByStructure groupAtomsByStructure(
    const std::vector<int>& structureTypes,
    const std::function<bool(std::size_t atomIndex)>& inCluster
);

// The per-structure atoms table: the counts in main_listing, the structures in sub_listings, and in
// export.AtomisticExporter the atoms of each structure under its export name, each with its identifier, position and
// structure, and the fields of the analysis atomFields adds.
json atomsByStructureTable(
    const LammpsParser::Frame& frame,
    const AtomsByStructure& structures,
    const std::function<void(json& atom, std::size_t atomIndex)>& atomFields
);

using StructureNameResolver = std::function<std::string(std::size_t atomIndex, int structureType)>;

using AtomColumnWriter = std::function<void(ColumnarAtomWriter& writer, std::size_t atomIndex, int structureType)>;

void streamStructureIdentificationToParquet(
    const std::string& filePath,
    const LammpsParser::Frame& frame,
    const StructureAnalysis& analysis,
    StructureNameResolver resolveStructureName = {},
    AtomColumnWriter atomColumnWriter = {}
);

// Writes the structure of each atom next to the output base, as <base>_atoms.parquet.
AtomsTable writeAtomsTable(
    const std::string& outputBase,
    const LammpsParser::Frame& frame,
    const StructureAnalysis& analysis,
    StructureNameResolver resolveStructureName = {},
    AtomColumnWriter atomColumnWriter = {}
);

}
