#include <volt/analysis/structure_identification_export.h>

#include <map>
#include <string>

#include <volt/structures/crystal_structure_types.h>

namespace Volt::StructureIdentificationExport {

namespace {

int atomIdForExport(const LammpsParser::Frame& frame, std::size_t atomIndex){
    return atomIndex < frame.ids.size()
        ? frame.ids[atomIndex]
        : static_cast<int>(atomIndex);
}

Point3 atomPositionForExport(const LammpsParser::Frame& frame, std::size_t atomIndex){
    if(atomIndex < frame.positions.size()){
        return frame.positions[atomIndex];
    }

    return Point3::Origin();
}

std::string defaultStructureName(int structureType){
    return structureTypeName(structureType);
}

std::string topologyNameForAtomExport(
    const StructureAnalysis& analysis,
    std::size_t atomIndex
){
    if(const Cluster* cluster = analysis.atomCluster(static_cast<int>(atomIndex));
       cluster && !cluster->topologyName.empty()){
        return cluster->topologyName;
    }

    return {};
}

json buildAtomRecord(
    const LammpsParser::Frame& frame,
    const StructureAnalysis& analysis,
    std::size_t atomIndex,
    int structureType,
    const std::string& structureName,
    const AtomRecordAugmenter& augmentAtom
){
    const StructureContext& context = analysis.context();
    const Point3 position = atomPositionForExport(frame, atomIndex);
    json atom = {
        {"id", atomIdForExport(frame, atomIndex)},
        {"pos", {position.x(), position.y(), position.z()}},
        {"structure_id", structureType},
        {"structure_name", structureName},
        {"cluster_id", context.atomClusters ? context.atomClusters->getInt(atomIndex) : 0}
    };

    const std::string topologyName = topologyNameForAtomExport(analysis, atomIndex);
    if(!topologyName.empty()){
        atom["topology_name"] = topologyName;
    }

    if(augmentAtom){
        augmentAtom(atom, atomIndex, structureType);
    }

    return atom;
}

struct StructureGroup{
    int structureId = static_cast<int>(StructureType::OTHER);
    std::string structureName;
    int atomCount = 0;
    json atoms = json::array();
};

}

json buildStructureIdentificationJson(
    const LammpsParser::Frame& frame,
    const StructureAnalysis& analysis,
    StructureNameResolver resolveStructureName,
    AtomRecordAugmenter augmentAtom
){
    const StructureContext& context = analysis.context();
    std::map<std::string, StructureGroup> groupsByName;
    json perAtomProperties = json::array();
    int clusteredAtoms = 0;

    for(std::size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
        const int structureType = context.structureTypes
            ? context.structureTypes->getInt(atomIndex)
            : static_cast<int>(StructureType::OTHER);
        const std::string structureName = resolveStructureName
            ? resolveStructureName(atomIndex, structureType)
            : defaultStructureName(structureType);
        auto& group = groupsByName[structureName];
        group.structureId = structureType;
        group.structureName = structureName;
        group.atomCount += 1;
        group.atoms.push_back(buildAtomRecord(
            frame,
            analysis,
            atomIndex,
            structureType,
            structureName,
            augmentAtom
        ));

        if(context.atomClusters && context.atomClusters->getInt(atomIndex) != 0){
            ++clusteredAtoms;
        }

        json atomProperties = {
            {"id", atomIdForExport(frame, atomIndex)},
            {"structure_id", structureType},
            {"structure_name", structureName},
            {"cluster_id", context.atomClusters ? context.atomClusters->getInt(atomIndex) : 0}
        };
        const std::string topologyName = topologyNameForAtomExport(analysis, atomIndex);
        if(!topologyName.empty()){
            atomProperties["topology_name"] = topologyName;
        }
        perAtomProperties.push_back(std::move(atomProperties));
    }

    json structures = json::array();
    json atomsByStructure = json::object();
    json chartNames = json::array();
    json chartCounts = json::array();
    for(auto& [structureName, group] : groupsByName){
        if(group.atomCount <= 0){
            continue;
        }

        structures.push_back({
            {"structure_id", group.structureId},
            {"structure_name", structureName},
            {"atom_count", group.atomCount}
        });
        chartNames.push_back(structureName);
        chartCounts.push_back(group.atomCount);
        atomsByStructure[structureName] = std::move(group.atoms);
    }

    json result;
    result["main_listing"] = {
        {"total_atoms", static_cast<int>(context.atomCount())},
        {"structure_count", static_cast<int>(structures.size())},
        {"clustered_atoms", clusteredAtoms},
        {"unclustered_atoms", static_cast<int>(context.atomCount()) - clusteredAtoms}
    };
    result["sub_listings"] = {
        {"structures", structures}
    };
    result["per-atom-properties"] = perAtomProperties;
    result["export"]["AtomisticExporter"] = std::move(atomsByStructure);
    result["export"]["ChartExporter"]["structure_counts"] = {
        {"structure_name", chartNames},
        {"atom_count", chartCounts}
    };
    return result;
}

}
