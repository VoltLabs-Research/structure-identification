#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/reconstructed_dump_utils.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/structures/crystal_structure_types.h>

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <duckdb.hpp>

namespace Volt{

namespace AnalysisContextDetail{

std::shared_ptr<ParticleProperty> makeNeighborIndicesProperty(const AnalysisContext& context){
    const std::size_t atomCount = context.atomCount();
    auto property = std::make_shared<ParticleProperty>(
        atomCount,
        DataType::Int,
        MAX_NEIGHBORS,
        0,
        true
    );
    std::fill(
        property->dataInt(),
        property->dataInt() + atomCount * MAX_NEIGHBORS,
        -1
    );

    if(!context.neighborCounts || !context.neighborOffsets || !context.neighborIndices){
        return property;
    }

    const int* counts = context.neighborCounts->constDataInt();
    const int* offsets = context.neighborOffsets->constDataInt();
    const int* indices = context.neighborIndices->constDataInt();
    const std::size_t totalNeighborEntries = context.neighborIndices->size();
    int* destination = property->dataInt();

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, atomCount), [&](const tbb::blocked_range<std::size_t>& range){
        for(std::size_t atomIndex = range.begin(); atomIndex < range.end(); ++atomIndex){
            const int count = std::max(0, std::min(counts[atomIndex], static_cast<int>(MAX_NEIGHBORS)));
            const int start = offsets[atomIndex];
            if(start < 0){
                continue;
            }
            for(int neighborSlot = 0; neighborSlot < count; ++neighborSlot){
                const std::size_t sourceIndex = static_cast<std::size_t>(start + neighborSlot);
                if(sourceIndex >= totalNeighborEntries){
                    break;
                }
                destination[atomIndex * MAX_NEIGHBORS + neighborSlot] = indices[sourceIndex];
            }
        }
    });

    return property;
}

std::shared_ptr<ParticleProperty> makeNeighborLatticeVectorProperty(
    const AnalysisContext& context,
    const StructureAnalysis& analysis
){
    auto property = std::make_shared<ParticleProperty>(
        context.atomCount(),
        DataType::Double,
        MAX_NEIGHBORS * 3,
        0.0,
        true
    );
    auto* data = property->dataDouble();
    const int* structureTypes = context.structureTypes->constDataInt();
    const int* neighborCounts = context.neighborCounts ? context.neighborCounts->constDataInt() : nullptr;

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, context.atomCount()), [&](const tbb::blocked_range<std::size_t>& range){
        for(std::size_t atomIndex = range.begin(); atomIndex < range.end(); ++atomIndex){
            const int structureType = structureTypes[atomIndex];
            if(structureType == LATTICE_OTHER){
                continue;
            }

            const int neighborCount = neighborCounts ? neighborCounts[atomIndex] : 0;
            const int maxCount = analysis.hasNeighborLatticeVectorOverrides()
                ? static_cast<int>(MAX_NEIGHBORS)
                : analysis.coordinationNumber(structureType);
            const int exportableCount = std::max(0, std::min(neighborCount, maxCount));
            const std::size_t atomBase = atomIndex * static_cast<std::size_t>(MAX_NEIGHBORS) * 3;

            for(int neighborSlot = 0; neighborSlot < exportableCount; ++neighborSlot){
                const Vector3& vector = analysis.neighborLatticeVector(static_cast<int>(atomIndex), neighborSlot);
                const std::size_t baseComponent = atomBase + static_cast<std::size_t>(neighborSlot) * 3;
                data[baseComponent + 0] = vector.x();
                data[baseComponent + 1] = vector.y();
                data[baseComponent + 2] = vector.z();
            }
        }
    });

    return property;
}

}

StructureContext::StructureContext(
    ParticleProperty* pos,
    const SimulationCell& cell,
    LatticeStructureType crystalType,
    ParticleProperty* outputStructures,
    std::vector<Matrix3>&& preferredOrientations
) : 
    positions(pos),
    structureTypes(outputStructures),
    maximumNeighborDistance(0.0),
    simCell(cell),
    inputCrystalType(crystalType),
    preferredCrystalOrientations(std::move(preferredOrientations))
{
    if(!positions){
        throw std::invalid_argument("Invalid positions");
    }
}

using namespace AnalysisContextDetail;

AnalysisContext::AnalysisContext(
    ParticleProperty* pos,
    const SimulationCell& cell,
    LatticeStructureType crystalType,
    ParticleProperty* outputStructures,
    std::vector<Matrix3>&& preferredOrientations
) :
    StructureContext(
        pos,
        cell,
        crystalType,
        outputStructures,
        std::move(preferredOrientations)
    )
{
    if(!structureTypes){
        throw std::invalid_argument("Invalid structure types");
    }
    const size_t numAtoms = atomCount();
    atomClusters = std::make_shared<ParticleProperty>(numAtoms, DataType::Int, 1, 0, true);
    atomAllowedSymmetryMasks = std::make_shared<ParticleProperty>(numAtoms, DataType::Int64, 1, 0, true);
    atomSymmetryPermutations = std::make_shared<ParticleProperty>(numAtoms, DataType::Int, 1, 0, true);
    std::fill(
        atomSymmetryPermutations->dataInt(),
        atomSymmetryPermutations->dataInt() + numAtoms,
        -1
    );

    if(numAtoms > 0){
        std::fill(
            structureTypes->dataInt(),
            structureTypes->dataInt() + numAtoms,
            LATTICE_OTHER
        );
    }
}

AnalysisContext::ExportedContext AnalysisContext::exportContext(
    const StructureAnalysis* analysis,
    const std::vector<ExtraScalarColumn>& extraColumns
) const{
    ExportedContext exported;
    auto& columns = exported.columns;
    auto& ownedProperties = exported.ownedProperties;
    const std::size_t count = atomCount();

    auto writeIntColumn = [&](const char* name, const std::shared_ptr<ParticleProperty>& property, int fillValue){
        auto exportedProperty = property ? property : AnalysisDumpUtils::makeIntProperty(count, fillValue);
        if(!property){
            ownedProperties.push_back(exportedProperty);
        }
        LammpsParser::writeColumn(columns, { name }, exportedProperty);
    };
    auto writeRawIntColumn = [&](const char* name, const ParticleProperty* property, int fillValue){
        if(!property){
            writeIntColumn(name, nullptr, fillValue);
            return;
        }

        auto exportedProperty = std::make_shared<ParticleProperty>(count, DataType::Int, 1, 0, true);
        std::copy(
            property->constDataInt(),
            property->constDataInt() + count,
            exportedProperty->dataInt()
        );
        ownedProperties.push_back(exportedProperty);
        LammpsParser::writeColumn(columns, { name }, exportedProperty);
    };
    writeRawIntColumn("structure_type", structureTypes, LATTICE_OTHER);
    writeIntColumn("cluster_id", atomClusters, 0);

    // Neighbor topology (neighbor_indices_* / neighbor_lattice_*) no longer rides
    // in the annotated dump — it is written to a dedicated `_neighbor_lattice.parquet`
    // sidecar by streamNeighborTopologyToParquet and consumed via inferFromContext.

    for(const auto& extra : extraColumns){
        if(!extra.property){
            continue;
        }
        ownedProperties.push_back(extra.property);
        LammpsParser::writeColumn(columns, { extra.name }, extra.property);
    }

    return exported;
}

bool AnalysisContext::writeDumpWithContext(
    const LammpsParser::Frame& frame,
    const std::string& outputFilename,
    const StructureAnalysis* analysis,
    const std::vector<ExtraScalarColumn>& extraColumns
) const{
    LammpsParser parser;
    const auto exported = exportContext(analysis, extraColumns);

    const std::vector<int>* atomIds = nullptr;
    std::vector<int> generatedAtomIds;

    if(frame.ids.size() == atomCount()){
        atomIds = &frame.ids;
    }else{
        generatedAtomIds.resize(atomCount());
        for(std::size_t atomIndex = 0; atomIndex < atomCount(); ++atomIndex){
            generatedAtomIds[atomIndex] = static_cast<int>(atomIndex);
        }
        atomIds = &generatedAtomIds;
    }

    return parser.writeFileMergedWithExtraColumns(
        outputFilename,
        frame,
        *atomIds,
        exported.columns,
        exported.headers,
        true
    );
}

namespace{

// Escapes a path for a single-quoted SQL string literal (mirrors CoreToolkit's
// internal Detail::sqlQuote, which plugins are not meant to include).
std::string sqlQuotePath(const std::string& path){
    std::string out;
    out.reserve(path.size() + 2);
    out.push_back('\'');
    for(char c : path){
        if(c == '\'') out.push_back('\'');
        out.push_back(c);
    }
    out.push_back('\'');
    return out;
}

}

bool streamNeighborTopologyToParquet(
    const std::string& filePath,
    const LammpsParser::Frame& frame,
    const AnalysisContext& context,
    const StructureAnalysis& analysis
){
    using namespace AnalysisContextDetail;

    const std::size_t atomCount = context.atomCount();

    // Build both packed properties once (indices: 18 components/atom,
    // lattice vectors: 54 = MAX_NEIGHBORS*3 components/atom).
    const auto neighborIndices = makeNeighborIndicesProperty(context);
    const auto neighborLattice = makeNeighborLatticeVectorProperty(context, analysis);
    const int* indexData = neighborIndices->constDataInt();
    const double* latticeData = neighborLattice->constDataDouble();

    const std::vector<int>* atomIds = nullptr;
    std::vector<int> generatedAtomIds;
    if(frame.ids.size() == atomCount){
        atomIds = &frame.ids;
    }else{
        generatedAtomIds.resize(atomCount);
        for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
            generatedAtomIds[atomIndex] = static_cast<int>(atomIndex);
        }
        atomIds = &generatedAtomIds;
    }

    try{
        // ponytail: fixed 74-col schema, Appender streams row-at-a-time — no
        // per-cell duckdb::Value buffer blowup like the dynamic atom writer.
        duckdb::DuckDB db(nullptr);
        duckdb::Connection con(db);

        std::string ddl = "CREATE TABLE neighbors(id UBIGINT, atom_index UINTEGER";
        for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
            ddl += ", \"" + AnalysisDumpUtils::neighborIndexName(slot) + "\" INTEGER";
        }
        for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
            ddl += ", \"" + AnalysisDumpUtils::neighborLatticeComponentName('x', slot) + "\" DOUBLE";
            ddl += ", \"" + AnalysisDumpUtils::neighborLatticeComponentName('y', slot) + "\" DOUBLE";
            ddl += ", \"" + AnalysisDumpUtils::neighborLatticeComponentName('z', slot) + "\" DOUBLE";
        }
        ddl += ')';
        if(con.Query(ddl)->HasError()){
            return false;
        }

        {
            duckdb::Appender appender(con, "neighbors");
            for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
                appender.BeginRow();
                appender.Append<std::uint64_t>(static_cast<std::uint64_t>((*atomIds)[atomIndex]));
                appender.Append<std::uint32_t>(static_cast<std::uint32_t>(atomIndex));

                const std::size_t indexBase = atomIndex * static_cast<std::size_t>(MAX_NEIGHBORS);
                for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
                    appender.Append<std::int32_t>(indexData[indexBase + static_cast<std::size_t>(slot)]);
                }

                const std::size_t latticeBase = atomIndex * static_cast<std::size_t>(MAX_NEIGHBORS) * 3;
                for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
                    const std::size_t component = latticeBase + static_cast<std::size_t>(slot) * 3;
                    appender.Append<double>(latticeData[component + 0]);
                    appender.Append<double>(latticeData[component + 1]);
                    appender.Append<double>(latticeData[component + 2]);
                }
                appender.EndRow();
            }
            appender.Close();
        }

        const std::string copySql =
            "COPY neighbors TO " + sqlQuotePath(filePath) +
            " (FORMAT PARQUET, COMPRESSION ZSTD)";
        return !con.Query(copySql)->HasError();
    }catch(const std::exception&){
        return false;
    }
}

}
