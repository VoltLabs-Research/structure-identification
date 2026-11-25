#include <volt/analysis/structure_identification_export.h>

#include <algorithm>
#include <string>
#include <string_view>
#include <utility>

#include <volt/plugin/output_files.h>
#include <volt/structures/crystal_structure_types.h>

namespace Volt::StructureIdentificationExport {

std::string defaultStructureName(int structureType){
    return structureTypeName(structureType);
}

const char* structureTypeNameForExport(int structureType){
    switch(static_cast<StructureType>(structureType)){
        case StructureType::SC:
            return "SC";
        case StructureType::FCC:
            return "FCC";
        case StructureType::HCP:
            return "HCP";
        case StructureType::BCC:
            return "BCC";
        case StructureType::CUBIC_DIAMOND:
            return "CUBIC_DIAMOND";
        case StructureType::HEX_DIAMOND:
            return "HEX_DIAMOND";
        case StructureType::ICO:
            return "ICO";
        case StructureType::GRAPHENE:
            return "GRAPHENE";
        case StructureType::CUBIC_DIAMOND_FIRST_NEIGH:
            return "CUBIC_DIAMOND_FIRST_NEIGH";
        case StructureType::CUBIC_DIAMOND_SECOND_NEIGH:
            return "CUBIC_DIAMOND_SECOND_NEIGH";
        case StructureType::HEX_DIAMOND_FIRST_NEIGH:
            return "HEX_DIAMOND_FIRST_NEIGH";
        case StructureType::HEX_DIAMOND_SECOND_NEIGH:
            return "HEX_DIAMOND_SECOND_NEIGH";
        case StructureType::OTHER:
        case StructureType::NUM_STRUCTURE_TYPES:
        default:
            return "OTHER";
    }
}

AtomsByStructure groupAtomsByStructure(
    const std::vector<int>& structureTypes,
    const std::function<bool(std::size_t atomIndex)>& inCluster
){
    constexpr int K = static_cast<int>(StructureType::NUM_STRUCTURE_TYPES);
    AtomsByStructure structures;
    structures.atoms = structureTypes.size();
    structures.atomsOfType.resize(K);
    for(std::size_t i = 0; i < structures.atoms; ++i){
        const int raw = structureTypes[i];
        const int st = (0 <= raw && raw < K) ? raw : 0;
        structures.atomsOfType[static_cast<std::size_t>(st)].push_back(i);
    }

    structures.order.reserve(K);
    for(int st = 0; st < K; st++){
        if(!structures.atomsOfType[static_cast<std::size_t>(st)].empty())
            structures.order.push_back(st);
    }
    std::sort(structures.order.begin(), structures.order.end(), [](int a, int b){
        return std::string_view(structureTypeNameForExport(a)) < std::string_view(structureTypeNameForExport(b));
    });

    for(std::size_t i = 0; i < structures.atoms; ++i){
        if(inCluster(i)) ++structures.clusteredAtoms;
    }
    return structures;
}

json atomsByStructureTable(
    const LammpsParser::Frame& frame,
    const AtomsByStructure& structures,
    const std::function<void(json& atom, std::size_t atomIndex)>& atomFields
){
    json structuresListing = json::array();
    json atomsByStructure;
    for(int st : structures.order){
        const char* name = structureTypeNameForExport(st);
        const std::vector<std::size_t>& atomsOfType = structures.atomsOfType[static_cast<std::size_t>(st)];
        json atomsArray = json::array();
        for(std::size_t atomIndex : atomsOfType){
            const Point3& pos = frame.positions[atomIndex];
            json atom = {
                {"id", frame.ids[atomIndex]},
                {"pos", {pos.x(), pos.y(), pos.z()}},
                {"structure_id", st},
                {"structure_name", name}
            };
            atomFields(atom, atomIndex);
            atomsArray.push_back(std::move(atom));
        }
        atomsByStructure[name] = std::move(atomsArray);
        structuresListing.push_back({
            {"structure_id", st},
            {"structure_name", name},
            {"atom_count", static_cast<int>(atomsOfType.size())}
        });
    }

    const int totalAtoms = frame.natoms;
    json table;
    table["main_listing"] = {
        {"total_atoms", totalAtoms},
        {"structure_count", static_cast<int>(structures.order.size())},
        {"clustered_atoms", structures.clusteredAtoms},
        {"unclustered_atoms", totalAtoms - structures.clusteredAtoms}
    };
    table["sub_listings"] = {
        {"structures", structuresListing}
    };
    table["export"] = json::object();
    table["export"]["AtomisticExporter"] = std::move(atomsByStructure);
    return table;
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

AtomsTable writeAtomsTable(
    const std::string& outputBase,
    const LammpsParser::Frame& frame,
    const StructureAnalysis& analysis,
    StructureNameResolver resolveStructureName,
    AtomColumnWriter atomColumnWriter
){
    AtomsTable atomsTable{outputPath(outputBase, "atoms")};
    streamStructureIdentificationToParquet(
        atomsTable.path, frame, analysis, std::move(resolveStructureName), std::move(atomColumnWriter)
    );
    return atomsTable;
}

}
