#include <volt/analysis/structure_identification_export.h>

#include <string>

#include <volt/structures/crystal_structure_types.h>

namespace Volt::StructureIdentificationExport {

std::string defaultStructureName(int structureType){
    return structureTypeName(structureType);
}

void streamStructureIdentificationToParquet(
    const std::string& filePath,
    const LammpsParser::Frame& frame,
    const StructureAnalysis& analysis,
    StructureNameResolver resolveStructureName,
    AtomColumnWriter atomColumnWriter
){
    const StructureContext& context = analysis.context();
    const int* structureTypes = context.structureTypes
        ? context.structureTypes->constDataInt()
        : nullptr;
    const int* atomClusters = context.atomClusters
        ? context.atomClusters->constDataInt()
        : nullptr;

    auto structureTypeFor = [&](std::size_t i) -> int {
        return structureTypes ? structureTypes[i] : static_cast<int>(StructureType::OTHER);
    };

    streamAtomsToParquet(
        filePath,
        frame,
        [&](std::size_t i){
            const int stype = structureTypeFor(i);
            return resolveStructureName ? resolveStructureName(i, stype) : defaultStructureName(stype);
        },
        [&](ColumnarAtomWriter& w, std::size_t i){
            w.field("cluster_id", atomClusters ? atomClusters[i] : 0);
            if(atomColumnWriter) atomColumnWriter(w, i, structureTypeFor(i));
        },
        [&](std::size_t i){ return structureTypeFor(i); },
        /*includeStructureColumns=*/true
    );
}

}
