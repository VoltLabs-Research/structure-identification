#pragma once

#include <cstddef>
#include <functional>
#include <string>

#include <volt/utilities/parquet_atom_writer.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/core/lammps_parser.h>

namespace Volt::StructureIdentificationExport {

using StructureNameResolver = std::function<std::string(std::size_t atomIndex, int structureType)>;

// Emits plugin-specific per-atom columns (e.g. PTM rmsd/orientation). The
// structure type for the atom is provided for convenience.
using AtomColumnWriter = std::function<void(ColumnarAtomWriter& writer, std::size_t atomIndex, int structureType)>;

// Writes the per-atom Parquet table for a structure-identification analysis:
// fixed columns (id, x/y/z, bucket=structure_name, structure_id=StructureType)
// plus cluster_id, optional topology_name, and any plugin columns.
// Listings/structure counts are derivable via GROUP BY structure_name.
void streamStructureIdentificationToParquet(
    const std::string& filePath,
    const LammpsParser::Frame& frame,
    const StructureAnalysis& analysis,
    StructureNameResolver resolveStructureName = {},
    AtomColumnWriter atomColumnWriter = {}
);

}
