#pragma once

#include <functional>
#include <string>

#include <nlohmann/json.hpp>

#include <volt/analysis/structure_analysis.h>
#include <volt/core/lammps_parser.h>

namespace Volt::StructureIdentificationExport {

using json = nlohmann::json;
using StructureNameResolver = std::function<std::string(std::size_t atomIndex, int structureType)>;
using AtomRecordAugmenter = std::function<void(json& atom, std::size_t atomIndex, int structureType)>;

json buildStructureIdentificationJson(
    const LammpsParser::Frame& frame,
    const StructureAnalysis& analysis,
    StructureNameResolver resolveStructureName = {},
    AtomRecordAugmenter augmentAtom = {}
);

}
