#pragma once

#include <cstddef>
#include <functional>
#include <string>

#include <volt/utilities/parquet_atom_writer.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/core/lammps_parser.h>

namespace Volt::StructureIdentificationExport {

using StructureNameResolver = std::function<std::string(std::size_t atomIndex, int structureType)>;

using AtomColumnWriter = std::function<void(ColumnarAtomWriter& writer, std::size_t atomIndex, int structureType)>;

void streamStructureIdentificationToParquet(
    const std::string& filePath,
    const LammpsParser::Frame& frame,
    const StructureAnalysis& analysis,
    StructureNameResolver resolveStructureName = {},
    AtomColumnWriter atomColumnWriter = {}
);

}
