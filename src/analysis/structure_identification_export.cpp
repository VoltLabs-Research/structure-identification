#include <volt/analysis/structure_identification_export.h>

#include <map>
#include <string>
#include <fstream>

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


}

void streamStructureIdentificationToFile(
    const std::string& filePath,
    const LammpsParser::Frame& frame,
    const StructureAnalysis& analysis,
    StructureNameResolver resolveStructureName,
    AtomFieldWriter atomFieldWriter
){
    const StructureContext& context = analysis.context();
    const std::size_t natoms = context.atomCount();

    // Pass 1: count atoms per structure and compute stats
    std::map<std::string, std::pair<int,int>> structureCounts; // name → (structureId, count)
    int clusteredAtoms = 0;
    for(std::size_t i = 0; i < natoms; ++i){
        const int stype = context.structureTypes
            ? context.structureTypes->getInt(i)
            : static_cast<int>(StructureType::OTHER);
        const std::string sname = resolveStructureName
            ? resolveStructureName(i, stype)
            : defaultStructureName(stype);
        structureCounts[sname].first = stype;
        structureCounts[sname].second++;
        if(context.atomClusters && context.atomClusters->getInt(i) != 0)
            ++clusteredAtoms;
    }

    // Determine extra fields per atom from the writer callback (use atom 0 as probe)
    int extraFields = 0;
    if(atomFieldWriter && natoms > 0){
        // Probe with a null sink to count fields
        struct NullBuf : std::streambuf{
            int overflow(int c) override { return c; }
        } nullBuf;
        std::ostream nullStream(&nullBuf);
        MsgpackWriter probe(nullStream);
        const int stype0 = context.structureTypes ? context.structureTypes->getInt(0) : 0;
        atomFieldWriter(probe, 0, stype0, extraFields);
    }
    // Base per-atom fields: id, pos, structure_id, structure_name, cluster_id
    const int baseFields = 5;

    std::ofstream of(filePath, std::ios::binary);
    MsgpackWriter w(of);

    // Root map: export, main_listing, sub_listings, per-atom-properties
    w.write_map_header(4);

    // "export"
    w.write_key("export");
    w.write_map_header(2);
    {
        // AtomisticExporter: atoms grouped by structure
        w.write_key("AtomisticExporter");
        w.write_map_header(static_cast<uint32_t>(structureCounts.size()));
        for(const auto& [sname, sc] : structureCounts){
            w.write_key(sname);
            w.write_array_header(static_cast<uint32_t>(sc.second));
            // Write atoms for this structure
            for(std::size_t i = 0; i < natoms; ++i){
                const int stype = context.structureTypes
                    ? context.structureTypes->getInt(i)
                    : static_cast<int>(StructureType::OTHER);
                const std::string iname = resolveStructureName
                    ? resolveStructureName(i, stype)
                    : defaultStructureName(stype);
                if(iname != sname) continue;

                int extra = extraFields;
                const std::string topo = topologyNameForAtomExport(analysis, i);
                if(!topo.empty()) extra++;

                w.write_map_header(static_cast<uint32_t>(baseFields + extra));
                w.write_key("id");
                w.write_int(i < frame.ids.size() ? frame.ids[i] : static_cast<int>(i));
                w.write_key("pos");
                w.write_array_header(3);
                const auto& pos = i < frame.positions.size() ? frame.positions[i] : Point3::Origin();
                w.write_double(pos.x()); w.write_double(pos.y()); w.write_double(pos.z());
                w.write_key("structure_id"); w.write_int(stype);
                w.write_key("structure_name"); w.write_str(iname);
                w.write_key("cluster_id");
                w.write_int(context.atomClusters ? context.atomClusters->getInt(i) : 0);
                if(!topo.empty()){ w.write_key("topology_name"); w.write_str(topo); }
                if(atomFieldWriter){
                    int dummy = 0;
                    atomFieldWriter(w, i, stype, dummy);
                }
            }
        }

        // ChartExporter
        w.write_key("ChartExporter");
        w.write_map_header(1);
        w.write_key("structure_counts");
        w.write_map_header(2);
        w.write_key("structure_name");
        w.write_array_header(static_cast<uint32_t>(structureCounts.size()));
        for(const auto& [sname, _] : structureCounts) w.write_str(sname);
        w.write_key("atom_count");
        w.write_array_header(static_cast<uint32_t>(structureCounts.size()));
        for(const auto& [_, sc] : structureCounts) w.write_int(sc.second);
    }

    // "main_listing"
    w.write_key("main_listing");
    w.write_map_header(4);
    w.write_key("total_atoms"); w.write_int(static_cast<int64_t>(natoms));
    w.write_key("structure_count"); w.write_int(static_cast<int64_t>(structureCounts.size()));
    w.write_key("clustered_atoms"); w.write_int(clusteredAtoms);
    w.write_key("unclustered_atoms"); w.write_int(static_cast<int64_t>(natoms) - clusteredAtoms);

    // "sub_listings"
    w.write_key("sub_listings");
    w.write_map_header(1);
    w.write_key("structures");
    w.write_array_header(static_cast<uint32_t>(structureCounts.size()));
    for(const auto& [sname, sc] : structureCounts){
        w.write_map_header(3);
        w.write_key("structure_id"); w.write_int(sc.first);
        w.write_key("structure_name"); w.write_str(sname);
        w.write_key("atom_count"); w.write_int(sc.second);
    }

    // "per-atom-properties" — flat array, one entry per atom
    w.write_key("per-atom-properties");
    w.write_array_header(static_cast<uint32_t>(natoms));
    for(std::size_t i = 0; i < natoms; ++i){
        const int stype = context.structureTypes
            ? context.structureTypes->getInt(i)
            : static_cast<int>(StructureType::OTHER);
        const std::string sname = resolveStructureName
            ? resolveStructureName(i, stype)
            : defaultStructureName(stype);
        const std::string topo = topologyNameForAtomExport(analysis, i);
        const int topoExtra = topo.empty() ? 0 : 1;
        w.write_map_header(static_cast<uint32_t>(4 + topoExtra));
        w.write_key("id");
        w.write_int(i < frame.ids.size() ? frame.ids[i] : static_cast<int>(i));
        w.write_key("structure_id"); w.write_int(stype);
        w.write_key("structure_name"); w.write_str(sname);
        w.write_key("cluster_id");
        w.write_int(context.atomClusters ? context.atomClusters->getInt(i) : 0);
        if(!topo.empty()){ w.write_key("topology_name"); w.write_str(topo); }
    }

    of.flush();
}

}
