#pragma once

#include <functional>
#include <string>
#include <fstream>

#include <nlohmann/json.hpp>
#include <volt/utilities/json_utils.h>

#include <volt/analysis/structure_analysis.h>
#include <volt/core/lammps_parser.h>

namespace Volt::StructureIdentificationExport {

using json = nlohmann::json;
using StructureNameResolver = std::function<std::string(std::size_t atomIndex, int structureType)>;

// Streaming augmenter: writes extra fields directly to MsgpackWriter.
// extraFieldCount must be the exact number of extra key/value pairs written.
using AtomFieldWriter = std::function<void(MsgpackWriter& w, std::size_t atomIndex, int structureType, int& extraFieldCount)>;

void streamStructureIdentificationToFile(
    const std::string& filePath,
    const LammpsParser::Frame& frame,
    const StructureAnalysis& analysis,
    StructureNameResolver resolveStructureName = {},
    AtomFieldWriter atomFieldWriter = {}
);

}
