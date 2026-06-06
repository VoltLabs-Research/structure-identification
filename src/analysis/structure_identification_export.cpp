#include <volt/analysis/structure_identification_export.h>

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <volt/structures/crystal_structure_types.h>

namespace Volt::StructureIdentificationExport {

namespace {

struct StructureGroup{
    int structureId = static_cast<int>(StructureType::OTHER);
    std::vector<std::uint32_t> atomIndices;
};

struct ExportFrameView{
    const int* ids = nullptr;
    const Point3* positions = nullptr;
};

constexpr std::size_t kMsgpackChunkSize = 16384;
constexpr int kAtomisticBaseFields = 5;
constexpr int kPerAtomBaseFields = 4;

std::string defaultStructureName(int structureType){
    return structureTypeName(structureType);
}

int atomIdForExport(const ExportFrameView& frameView, std::size_t atomIndex){
    return frameView.ids
        ? frameView.ids[atomIndex]
        : static_cast<int>(atomIndex);
}

Point3 atomPositionForExport(const ExportFrameView& frameView, std::size_t atomIndex){
    return frameView.positions
        ? frameView.positions[atomIndex]
        : Point3::Origin();
}

void renderAtomisticChunk(
    std::string& renderedChunk,
    const StructureGroup& group,
    std::size_t chunkBegin,
    std::size_t chunkEnd,
    const ExportFrameView& frameView,
    const int* structureTypes,
    const int* atomClusters,
    const std::vector<std::string_view>& topologyNames,
    const std::string& structureName,
    const AtomFieldWriter& atomFieldWriter,
    const AtomFieldCountResolver& atomFieldCountResolver,
    int fixedExtraFields
){
    std::ostringstream stream(std::ios::binary | std::ios::out);
    MsgpackWriter writer(stream);

    for(std::size_t groupOffset = chunkBegin; groupOffset < chunkEnd; ++groupOffset){
        const std::size_t atomIndex = static_cast<std::size_t>(group.atomIndices[groupOffset]);
        const int structureType = structureTypes
            ? structureTypes[atomIndex]
            : static_cast<int>(StructureType::OTHER);
        const int clusterId = atomClusters ? atomClusters[atomIndex] : 0;
        int extraFields = atomFieldCountResolver
            ? atomFieldCountResolver(atomIndex, structureType)
            : fixedExtraFields;
        const std::string_view topologyName = topologyNames[atomIndex];
        if(!topologyName.empty()){
            ++extraFields;
        }

        writer.write_map_header(static_cast<uint32_t>(kAtomisticBaseFields + extraFields));
        writer.write_key("id");
        writer.write_int(atomIdForExport(frameView, atomIndex));
        writer.write_key("pos");
        writer.write_array_header(3);
        const Point3 position = atomPositionForExport(frameView, atomIndex);
        writer.write_double(position.x());
        writer.write_double(position.y());
        writer.write_double(position.z());
        writer.write_key("structure_id");
        writer.write_int(structureType);
        writer.write_key("structure_name");
        writer.write_str(structureName);
        writer.write_key("cluster_id");
        writer.write_int(clusterId);
        if(!topologyName.empty()){
            writer.write_key("topology_name");
            writer.write_str(topologyName);
        }
        if(atomFieldWriter){
            int writtenExtraFields = 0;
            atomFieldWriter(writer, atomIndex, structureType, writtenExtraFields);
        }
    }

    renderedChunk = stream.str();
}

void renderPerAtomChunk(
    std::string& renderedChunk,
    std::size_t atomBegin,
    std::size_t atomEnd,
    const ExportFrameView& frameView,
    const int* structureTypes,
    const int* atomClusters,
    const std::vector<std::uint32_t>& atomGroupOrdinals,
    const std::vector<const std::string*>& orderedStructureNames,
    const std::vector<std::string_view>& topologyNames
){
    std::ostringstream stream(std::ios::binary | std::ios::out);
    MsgpackWriter writer(stream);

    for(std::size_t atomIndex = atomBegin; atomIndex < atomEnd; ++atomIndex){
        const int structureType = structureTypes
            ? structureTypes[atomIndex]
            : static_cast<int>(StructureType::OTHER);
        const std::string& structureName = *orderedStructureNames[atomGroupOrdinals[atomIndex]];
        const std::string_view topologyName = topologyNames[atomIndex];
        const int extraFields = topologyName.empty() ? 0 : 1;

        writer.write_map_header(static_cast<uint32_t>(kPerAtomBaseFields + extraFields));
        writer.write_key("id");
        writer.write_int(atomIdForExport(frameView, atomIndex));
        writer.write_key("structure_id");
        writer.write_int(structureType);
        writer.write_key("structure_name");
        writer.write_str(structureName);
        writer.write_key("cluster_id");
        writer.write_int(atomClusters ? atomClusters[atomIndex] : 0);
        if(!topologyName.empty()){
            writer.write_key("topology_name");
            writer.write_str(topologyName);
        }
    }

    renderedChunk = stream.str();
}

}

void streamStructureIdentificationToFile(
    const std::string& filePath,
    const LammpsParser::Frame& frame,
    const StructureAnalysis& analysis,
    StructureNameResolver resolveStructureName,
    AtomFieldWriter atomFieldWriter,
    AtomFieldCountResolver atomFieldCountResolver
){
    const StructureContext& context = analysis.context();
    const std::size_t natoms = context.atomCount();
    const int* structureTypes = context.structureTypes
        ? context.structureTypes->constDataInt()
        : nullptr;
    const int* atomClusters = context.atomClusters
        ? context.atomClusters->constDataInt()
        : nullptr;
    ExportFrameView frameView;
    frameView.ids = frame.ids.size() == natoms ? frame.ids.data() : nullptr;
    frameView.positions = frame.positions.size() == natoms ? frame.positions.data() : nullptr;

    std::map<std::string, StructureGroup> structureGroups;
    int clusteredAtoms = 0;
    for(std::size_t i = 0; i < natoms; ++i){
        const int stype = structureTypes
            ? structureTypes[i]
            : static_cast<int>(StructureType::OTHER);
        const std::string sname = resolveStructureName
            ? resolveStructureName(i, stype)
            : defaultStructureName(stype);

        StructureGroup& group = structureGroups[sname];
        group.structureId = stype;
        group.atomIndices.push_back(static_cast<std::uint32_t>(i));
        if(atomClusters && atomClusters[i] != 0){
            ++clusteredAtoms;
        }
    }

    std::vector<const std::string*> orderedStructureNames;
    std::vector<const StructureGroup*> orderedStructureGroups;
    orderedStructureNames.reserve(structureGroups.size());
    orderedStructureGroups.reserve(structureGroups.size());

    std::vector<std::uint32_t> atomGroupOrdinals(natoms, 0);
    std::uint32_t groupOrdinal = 0;
    for(const auto& [name, group] : structureGroups){
        orderedStructureNames.push_back(&name);
        orderedStructureGroups.push_back(&group);
        for(const std::uint32_t atomIndex : group.atomIndices){
            atomGroupOrdinals[static_cast<std::size_t>(atomIndex)] = groupOrdinal;
        }
        ++groupOrdinal;
    }

    std::vector<std::string_view> topologyNames(natoms);
    for(std::size_t atomIndex = 0; atomIndex < natoms; ++atomIndex){
        if(const Cluster* cluster = analysis.atomCluster(static_cast<int>(atomIndex));
           cluster && !cluster->topologyName.empty()){
            topologyNames[atomIndex] = cluster->topologyName;
        }
    }

    int fixedExtraFields = 0;
    if(!atomFieldCountResolver && atomFieldWriter && natoms > 0){
        // Probe with a null sink to count fields
        struct NullBuf : std::streambuf{
            int overflow(int c) override { return c; }
        } nullBuf;
        std::ostream nullStream(&nullBuf);
        MsgpackWriter probe(nullStream);
        const int stype0 = structureTypes ? structureTypes[0] : 0;
        atomFieldWriter(probe, 0, stype0, fixedExtraFields);
    }

    std::ofstream of(filePath, std::ios::binary);
    std::vector<char> fileBuffer(1 << 20);
    of.rdbuf()->pubsetbuf(fileBuffer.data(), static_cast<std::streamsize>(fileBuffer.size()));
    MsgpackWriter w(of);
    const std::size_t maxWaveChunks = std::max<std::size_t>(1, std::thread::hardware_concurrency());

    // Root map: export, main_listing, sub_listings, per-atom-properties
    w.write_map_header(4);

    // "export"
    w.write_key("export");
    w.write_map_header(2);
    {
        // AtomisticExporter: atoms grouped by structure
        w.write_key("AtomisticExporter");
        w.write_map_header(static_cast<uint32_t>(orderedStructureGroups.size()));
        for(std::size_t groupIndex = 0; groupIndex < orderedStructureGroups.size(); ++groupIndex){
            const std::string& sname = *orderedStructureNames[groupIndex];
            const StructureGroup& group = *orderedStructureGroups[groupIndex];
            w.write_key(sname);
            w.write_array_header(static_cast<uint32_t>(group.atomIndices.size()));
            const std::size_t groupSize = group.atomIndices.size();
            for(std::size_t waveStart = 0; waveStart < groupSize; waveStart += kMsgpackChunkSize * maxWaveChunks){
                const std::size_t remainingAtoms = groupSize - waveStart;
                const std::size_t waveChunkCount = std::min(
                    maxWaveChunks,
                    (remainingAtoms + kMsgpackChunkSize - 1) / kMsgpackChunkSize
                );
                std::vector<std::string> renderedChunks(waveChunkCount);

                tbb::parallel_for(
                    tbb::blocked_range<std::size_t>(0, waveChunkCount),
                    [&](const tbb::blocked_range<std::size_t>& range){
                        for(std::size_t chunkIndex = range.begin(); chunkIndex < range.end(); ++chunkIndex){
                            const std::size_t chunkBegin = waveStart + chunkIndex * kMsgpackChunkSize;
                            const std::size_t chunkEnd = std::min(chunkBegin + kMsgpackChunkSize, groupSize);
                            renderAtomisticChunk(
                                renderedChunks[chunkIndex],
                                group,
                                chunkBegin,
                                chunkEnd,
                                frameView,
                                structureTypes,
                                atomClusters,
                                topologyNames,
                                sname,
                                atomFieldWriter,
                                atomFieldCountResolver,
                                fixedExtraFields
                            );
                        }
                    }
                );

                for(const std::string& chunk : renderedChunks){
                    of.write(chunk.data(), static_cast<std::streamsize>(chunk.size()));
                }
            }
        }

        // ChartExporter
        w.write_key("ChartExporter");
        w.write_map_header(1);
        w.write_key("structure_counts");
        w.write_map_header(2);
        w.write_key("structure_name");
        w.write_array_header(static_cast<uint32_t>(orderedStructureGroups.size()));
        for(const std::string* sname : orderedStructureNames){
            w.write_str(*sname);
        }
        w.write_key("atom_count");
        w.write_array_header(static_cast<uint32_t>(orderedStructureGroups.size()));
        for(const StructureGroup* group : orderedStructureGroups){
            w.write_int(static_cast<int>(group->atomIndices.size()));
        }
    }

    // "main_listing"
    w.write_key("main_listing");
    w.write_map_header(4);
    w.write_key("total_atoms"); w.write_int(static_cast<int64_t>(natoms));
    w.write_key("structure_count"); w.write_int(static_cast<int64_t>(orderedStructureGroups.size()));
    w.write_key("clustered_atoms"); w.write_int(clusteredAtoms);
    w.write_key("unclustered_atoms"); w.write_int(static_cast<int64_t>(natoms) - clusteredAtoms);

    // "sub_listings"
    w.write_key("sub_listings");
    w.write_map_header(1);
    w.write_key("structures");
    w.write_array_header(static_cast<uint32_t>(orderedStructureGroups.size()));
    for(std::size_t groupIndex = 0; groupIndex < orderedStructureGroups.size(); ++groupIndex){
        const std::string& sname = *orderedStructureNames[groupIndex];
        const StructureGroup& group = *orderedStructureGroups[groupIndex];
        w.write_map_header(3);
        w.write_key("structure_id"); w.write_int(group.structureId);
        w.write_key("structure_name"); w.write_str(sname);
        w.write_key("atom_count"); w.write_int(static_cast<int>(group.atomIndices.size()));
    }

    // "per-atom-properties" — flat array, one entry per atom
    w.write_key("per-atom-properties");
    w.write_array_header(static_cast<uint32_t>(natoms));
    for(std::size_t waveStart = 0; waveStart < natoms; waveStart += kMsgpackChunkSize * maxWaveChunks){
        const std::size_t remainingAtoms = natoms - waveStart;
        const std::size_t waveChunkCount = std::min(
            maxWaveChunks,
            (remainingAtoms + kMsgpackChunkSize - 1) / kMsgpackChunkSize
        );
        std::vector<std::string> renderedChunks(waveChunkCount);

        tbb::parallel_for(
            tbb::blocked_range<std::size_t>(0, waveChunkCount),
            [&](const tbb::blocked_range<std::size_t>& range){
                for(std::size_t chunkIndex = range.begin(); chunkIndex < range.end(); ++chunkIndex){
                    const std::size_t chunkBegin = waveStart + chunkIndex * kMsgpackChunkSize;
                    const std::size_t chunkEnd = std::min(chunkBegin + kMsgpackChunkSize, natoms);
                    renderPerAtomChunk(
                        renderedChunks[chunkIndex],
                        chunkBegin,
                        chunkEnd,
                        frameView,
                        structureTypes,
                        atomClusters,
                        atomGroupOrdinals,
                        orderedStructureNames,
                        topologyNames
                    );
                }
            }
        );

        for(const std::string& chunk : renderedChunks){
            of.write(chunk.data(), static_cast<std::streamsize>(chunk.size()));
        }
    }

    of.flush();
}

}
