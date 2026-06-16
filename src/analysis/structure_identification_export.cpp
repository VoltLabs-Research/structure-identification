#include <volt/analysis/structure_identification_export.h>

#include <string>
#include <string_view>
#include <vector>

#include <volt/structures/crystal_structure_types.h>

namespace Volt::StructureIdentificationExport {

namespace {

std::string defaultStructureName(int structureType){
    return structureTypeName(structureType);
}

}

void streamStructureIdentificationToParquet(
    const std::string& filePath,
    const LammpsParser::Frame& frame,
    const StructureAnalysis& analysis,
    StructureNameResolver resolveStructureName,
    AtomColumnWriter atomColumnWriter
){
    const StructureContext& context = analysis.context();
    const std::size_t natoms = context.atomCount();
    const int* structureTypes = context.structureTypes
        ? context.structureTypes->constDataInt()
        : nullptr;
    const int* atomClusters = context.atomClusters
        ? context.atomClusters->constDataInt()
        : nullptr;

    auto structureTypeFor = [&](std::size_t i) -> int {
        return structureTypes ? structureTypes[i] : static_cast<int>(StructureType::OTHER);
    };

    // Topology name per atom (sparse); emitted as a column when present.
    std::vector<std::string_view> topologyNames(natoms);
    for(std::size_t atomIndex = 0; atomIndex < natoms; ++atomIndex){
        if(const Cluster* cluster = analysis.atomCluster(static_cast<int>(atomIndex));
           cluster && !cluster->topologyName.empty()){
            topologyNames[atomIndex] = cluster->topologyName;
        }
    }

    streamAtomsToParquet(
        filePath,
        frame,
        [&](std::size_t i){
            const int stype = structureTypeFor(i);
            return resolveStructureName ? resolveStructureName(i, stype) : defaultStructureName(stype);
        },
        [&](ColumnarAtomWriter& w, std::size_t i){
            w.field("cluster_id", atomClusters ? atomClusters[i] : 0);
            if(i < topologyNames.size() && !topologyNames[i].empty()){
                w.field("topology_name", std::string(topologyNames[i]));
            }
            if(atomColumnWriter) atomColumnWriter(w, i, structureTypeFor(i));
        },
        [&](std::size_t i){ return structureTypeFor(i); },
        // This helper is the structural-identification path (PTM/ACNA/opendxa/
        // structure-identification): always emit structure_id/structure_name.
        /*includeStructureColumns=*/true
    );
}

}
