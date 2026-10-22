#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/reconstructed_dump_utils.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/core/progress.h>
#include <volt/structures/crystal_structure_types.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>

#include <spdlog/spdlog.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <duckdb.hpp>
#include <volt/utilities/duckdb_parquet.h>

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

    for(const auto& extra : extraColumns){
        if(!extra.property){
            continue;
        }
        ownedProperties.push_back(extra.property);
        LammpsParser::writeColumn(columns, { extra.name }, extra.property);
    }

    return exported;
}

int contractAtomId(const LammpsParser::Frame& frame, std::size_t atomIndex){
    return atomIndex < frame.ids.size()
        ? frame.ids[atomIndex]
        : static_cast<int>(atomIndex + 1);
}

LammpsParser::Frame maskFrame(const LammpsParser::Frame& frame, const std::vector<bool>& atomMask){
    const std::size_t sourceCount = frame.positions.size();
    const std::size_t keptCount = static_cast<std::size_t>(
        std::count(atomMask.begin(), atomMask.end(), true));

    LammpsParser::Frame subset;
    subset.timestep = frame.timestep;
    subset.natoms = static_cast<int>(keptCount);
    subset.simulationCell = frame.simulationCell;
    subset.headerOrder = frame.headerOrder;
    subset.headerProperties = frame.headerProperties;
    subset.atomColumnOrder = frame.atomColumnOrder;
    subset.atomColumnsScaled = frame.atomColumnsScaled;

    const bool hasTypes = frame.types.size() == sourceCount;
    const bool hasImages = frame.hasImageFlags() && frame.imageX.size() == sourceCount;
    subset.positions.reserve(keptCount);
    subset.ids.reserve(keptCount);
    if(hasTypes){
        subset.types.reserve(keptCount);
    }
    for(std::size_t atomIndex = 0; atomIndex < sourceCount; ++atomIndex){
        if(atomIndex >= atomMask.size() || !atomMask[atomIndex]){
            continue;
        }
        subset.positions.push_back(frame.positions[atomIndex]);
        subset.ids.push_back(contractAtomId(frame, atomIndex));
        if(hasTypes){
            subset.types.push_back(frame.types[atomIndex]);
        }
        if(hasImages){
            subset.imageX.push_back(frame.imageX[atomIndex]);
            subset.imageY.push_back(frame.imageY[atomIndex]);
            subset.imageZ.push_back(frame.imageZ[atomIndex]);
        }
    }

    for(const auto& [name, column] : frame.atomProperties){
        if(column.size() != sourceCount){
            continue;
        }
        LammpsParser::AtomColumn filtered;
        filtered.dataType = column.dataType;
        for(std::size_t atomIndex = 0; atomIndex < sourceCount; ++atomIndex){
            if(atomIndex >= atomMask.size() || !atomMask[atomIndex]){
                continue;
            }
            switch(column.dataType){
                case DataType::Int:
                    filtered.ints.push_back(column.ints[atomIndex]);
                    break;
                case DataType::Int64:
                    filtered.int64s.push_back(column.int64s[atomIndex]);
                    break;
                case DataType::Double:
                    filtered.doubles.push_back(column.doubles[atomIndex]);
                    break;
                case DataType::Void:
                    break;
            }
        }
        subset.atomProperties.emplace(name, std::move(filtered));
    }

    return subset;
}

bool AnalysisContext::writeDumpWithContext(
    const LammpsParser::Frame& frame,
    const std::string& outputFilename,
    const StructureAnalysis* analysis,
    const std::vector<ExtraScalarColumn>& extraColumns,
    const std::vector<bool>* atomMask
) const{
    const Stopwatch timer;
    LammpsParser parser;
    const auto exported = exportContext(analysis, extraColumns);

    std::vector<int> propertyAtomIds(atomCount());
    for(std::size_t atomIndex = 0; atomIndex < atomCount(); ++atomIndex){
        propertyAtomIds[atomIndex] = contractAtomId(frame, atomIndex);
    }

    const auto reportWritten = [&](bool written){
        if(written){
            spdlog::info("  Annotated dump written to {} in {}", outputFilename, timer.elapsed());
        }
        return written;
    };

    if(!atomMask){
        return reportWritten(parser.writeFileMergedWithExtraColumns(
            outputFilename,
            frame,
            propertyAtomIds,
            exported.columns,
            exported.headers,
            true
        ));
    }

    const LammpsParser::Frame subset = maskFrame(frame, *atomMask);
    return reportWritten(parser.writeFileMergedWithExtraColumns(
        outputFilename,
        subset,
        propertyAtomIds,
        exported.columns,
        exported.headers,
        true
    ));
}

bool streamNeighborTopologyToParquet(
    const std::string& filePath,
    const LammpsParser::Frame& frame,
    const AnalysisContext& context,
    const StructureAnalysis& analysis,
    const std::vector<bool>* atomMask
){
    using namespace AnalysisContextDetail;

    const Stopwatch timer;
    const std::size_t atomCount = context.atomCount();

    const auto neighborIndices = makeNeighborIndicesProperty(context);
    const auto neighborLattice = makeNeighborLatticeVectorProperty(context, analysis);
    const int* indexData = neighborIndices->constDataInt();
    const double* latticeData = neighborLattice->constDataDouble();

    std::vector<int> compactIndex(atomCount, -1);
    int nextCompactIndex = 0;
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        if(atomMask && (atomIndex >= atomMask->size() || !(*atomMask)[atomIndex])){
            continue;
        }
        compactIndex[atomIndex] = nextCompactIndex++;
    }

    try{
        auto db = Volt::Detail::openInMemoryDb();
        duckdb::Connection con(*db);

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
        auto ddlResult = con.Query(ddl);
        if(ddlResult->HasError()){
            spdlog::error("  Neighbor lattice export to {} failed: {}", filePath, ddlResult->GetError());
            return false;
        }

        {
            duckdb::Appender appender(con, "neighbors");
            ProgressReporter progress("  Writing neighbor lattice", atomCount);
            for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
                progress.update(atomIndex);
                if(compactIndex[atomIndex] < 0){
                    continue;
                }
                appender.BeginRow();
                appender.Append<std::uint64_t>(static_cast<std::uint64_t>(contractAtomId(frame, atomIndex)));
                appender.Append<std::uint32_t>(static_cast<std::uint32_t>(compactIndex[atomIndex]));

                const std::size_t indexBase = atomIndex * static_cast<std::size_t>(MAX_NEIGHBORS);
                std::array<int, MAX_NEIGHBORS> remapped{};
                for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
                    const int neighborIndex = indexData[indexBase + static_cast<std::size_t>(slot)];
                    remapped[static_cast<std::size_t>(slot)] = neighborIndex >= 0
                        ? compactIndex[static_cast<std::size_t>(neighborIndex)]
                        : -1;
                    appender.Append<std::int32_t>(remapped[static_cast<std::size_t>(slot)]);
                }

                const std::size_t latticeBase = atomIndex * static_cast<std::size_t>(MAX_NEIGHBORS) * 3;
                for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
                    const std::size_t component = latticeBase + static_cast<std::size_t>(slot) * 3;
                    const bool kept = remapped[static_cast<std::size_t>(slot)] >= 0;
                    appender.Append<double>(kept ? latticeData[component + 0] : 0.0);
                    appender.Append<double>(kept ? latticeData[component + 1] : 0.0);
                    appender.Append<double>(kept ? latticeData[component + 2] : 0.0);
                }
                appender.EndRow();
            }
            appender.Close();
        }

        const std::string copySql =
            "COPY neighbors TO " + Detail::sqlQuote(filePath) +
            " (FORMAT PARQUET, COMPRESSION ZSTD)";
        auto copyResult = con.Query(copySql);
        if(copyResult->HasError()){
            spdlog::error("  Neighbor lattice export to {} failed: {}", filePath, copyResult->GetError());
            return false;
        }
        spdlog::info("  Neighbor lattice written to {} ({} atoms) in {}", filePath, nextCompactIndex, timer.elapsed());
        return true;
    }catch(const std::exception& error){
        spdlog::error("  Neighbor lattice export to {} failed: {}", filePath, error.what());
        return false;
    }
}

}
