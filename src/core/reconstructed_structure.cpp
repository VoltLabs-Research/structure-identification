#include <volt/core/reconstructed_structure.h>

#include <chrono>
#include <volt/analysis/reconstructed_dump_utils.h>

#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include <duckdb.hpp>
#include <volt/utilities/duckdb_parquet.h>

namespace Volt{

namespace{

bool assembleNeighborGraph(
    const LammpsParser::Frame& frame,
    StructureAnalysis& analysis,
    ReconstructedStructureContext& context,
    const std::array<const int*, MAX_NEIGHBORS>& neighborSlotPtrs,
    const std::array<std::array<const double*, MAX_NEIGHBORS>, 3>& latticePtrs
){
    const size_t natoms = static_cast<size_t>(frame.natoms);

    std::vector<int> neighborCounts(natoms, 0);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, natoms, 8192),
        [&](const tbb::blocked_range<size_t>& r){
        for(size_t atomIndex = r.begin(); atomIndex < r.end(); ++atomIndex){
            int count = 0;
            for(int s = 0; s < MAX_NEIGHBORS; ++s){
                if(neighborSlotPtrs[s][atomIndex] < 0) break;
                ++count;
            }
            neighborCounts[atomIndex] = count;
        }
    });

    std::vector<int> compactOffsets(natoms + 1, 0);
    int totalNeighborEntries = 0;
    for(size_t i = 0; i < natoms; ++i){
        compactOffsets[i] = totalNeighborEntries;
        totalNeighborEntries += neighborCounts[i];
    }
    compactOffsets[natoms] = totalNeighborEntries;

    context.neighborCounts = AnalysisDumpUtils::makeIntProperty(neighborCounts);
    context.neighborOffsets = AnalysisDumpUtils::makeIntProperty(compactOffsets);

    auto compactNeighborIndices = std::make_shared<ParticleProperty>(
        static_cast<size_t>(totalNeighborEntries), DataType::Int, 1, 0, true
    );
    int* compactData = compactNeighborIndices->dataInt();

    tbb::parallel_for(tbb::blocked_range<size_t>(0, natoms, 8192),
        [&](const tbb::blocked_range<size_t>& r){
        for(size_t atomIndex = r.begin(); atomIndex < r.end(); ++atomIndex){
            const int count = neighborCounts[atomIndex];
            const int start = compactOffsets[atomIndex];
            for(int s = 0; s < count; ++s){
                compactData[start + s] = neighborSlotPtrs[s][atomIndex];
            }
        }
    });

    context.neighborIndices = compactNeighborIndices;

    context.maximumNeighborDistance = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, natoms, 8192),
        0.0,
        [&](const tbb::blocked_range<size_t>& r, double localMax) -> double{
            for(size_t atomIndex = r.begin(); atomIndex < r.end(); ++atomIndex){
                const int count = neighborCounts[atomIndex];
                const int start = compactOffsets[atomIndex];
                for(int s = 0; s < count; ++s){
                    const int neighbor = compactData[start + s];
                    const Vector3 delta = context.simCell.wrapVector(
                        frame.positions[static_cast<size_t>(neighbor)] -
                        frame.positions[atomIndex]
                    );
                    localMax = std::max(localMax, delta.length());
                }
            }
            return localMax;
        },
        [](double a, double b){ return std::max(a, b); }
    );

    std::vector<Vector3> neighborVectorOverrides(
        context.atomCount() * static_cast<size_t>(MAX_NEIGHBORS),
        Vector3::Zero()
    );

    tbb::parallel_for(tbb::blocked_range<size_t>(0, context.atomCount(), 8192),
        [&](const tbb::blocked_range<size_t>& r){
        for(size_t atomIndex = r.begin(); atomIndex < r.end(); ++atomIndex){
            for(int s = 0; s < MAX_NEIGHBORS; ++s){
                neighborVectorOverrides[atomIndex * static_cast<size_t>(MAX_NEIGHBORS) + static_cast<size_t>(s)] = Vector3(
                    latticePtrs[0][s][atomIndex],
                    latticePtrs[1][s][atomIndex],
                    latticePtrs[2][s][atomIndex]
                );
            }
        }
    });

    analysis.setNeighborLatticeVectorOverrides(
        std::move(neighborVectorOverrides),
        static_cast<size_t>(MAX_NEIGHBORS)
    );

    return true;
}

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

ReconstructedStructureContext::ReconstructedStructureContext(
    ParticleProperty* positions,
    const SimulationCell& cell
) :
    StructureContext(positions, cell){}

bool ReconstructedStructureContext::loadStructureAndClusterFromFrame(
    const LammpsParser::Frame& frame,
    ReconstructedStructureContext& context,
    std::string* errorMessage
){
    if(frame.findAtomProperty("structure_type")){
        std::vector<int> structureTypes;
        if(!AnalysisDumpUtils::extractIntegralColumn(frame, "structure_type", structureTypes, errorMessage)){
            return false;
        }
        context._ownedStructureTypes = AnalysisDumpUtils::makeIntProperty(structureTypes);
        context.structureTypes = context._ownedStructureTypes.get();
    }

    std::vector<int> clusterIds;
    if(!AnalysisDumpUtils::extractIntegralColumn(frame, "cluster_id", clusterIds, errorMessage)){
        return false;
    }
    context.atomClusters = AnalysisDumpUtils::makeIntProperty(clusterIds);

    return true;
}

namespace{

template<typename Dest, typename Convert>
bool readColumnInto(duckdb::Vector& vec, duckdb::idx_t count, const std::vector<int>& targetRow,
                    duckdb::idx_t chunkOffset, Dest* dest, Convert convert){
    duckdb::UnifiedVectorFormat fmt;
    vec.ToUnifiedFormat(count, fmt);
    for(duckdb::idx_t r = 0; r < count; ++r){
        const int target = targetRow[chunkOffset + r];
        if(target < 0) continue;
        const auto idx = fmt.sel->get_index(r);
        if(!fmt.validity.RowIsValid(idx)) continue;
        convert(fmt, idx, dest[target]);
    }
    return true;
}

template<typename Dest>
bool readIntegerColumn(duckdb::Vector& vec, duckdb::idx_t count, const std::vector<int>& targetRow,
                       duckdb::idx_t chunkOffset, Dest* dest, std::string* what){
    const auto assign = [&](auto tag) {
        using Physical = decltype(tag);
        return readColumnInto(vec, count, targetRow, chunkOffset, dest,
            [](const duckdb::UnifiedVectorFormat& f, duckdb::idx_t i, Dest& out){
                out = static_cast<Dest>(duckdb::UnifiedVectorFormat::GetData<Physical>(f)[i]);
            });
    };

    switch(vec.GetType().id()){
        case duckdb::LogicalTypeId::TINYINT:   return assign(int8_t{});
        case duckdb::LogicalTypeId::SMALLINT:  return assign(int16_t{});
        case duckdb::LogicalTypeId::INTEGER:   return assign(int32_t{});
        case duckdb::LogicalTypeId::BIGINT:    return assign(int64_t{});
        case duckdb::LogicalTypeId::UTINYINT:  return assign(uint8_t{});
        case duckdb::LogicalTypeId::USMALLINT: return assign(uint16_t{});
        case duckdb::LogicalTypeId::UINTEGER:  return assign(uint32_t{});
        case duckdb::LogicalTypeId::UBIGINT:   return assign(uint64_t{});
        default:
            if(what) *what = "expected an integer column, got " + vec.GetType().ToString();
            return false;
    }
}

bool readDoubleColumn(duckdb::Vector& vec, duckdb::idx_t count, const std::vector<int>& targetRow,
                      duckdb::idx_t chunkOffset, double* dest, std::string* what){
    switch(vec.GetType().id()){
        case duckdb::LogicalTypeId::DOUBLE:
            return readColumnInto(vec, count, targetRow, chunkOffset, dest,
                [](const duckdb::UnifiedVectorFormat& f, duckdb::idx_t i, double& out){
                    out = duckdb::UnifiedVectorFormat::GetData<double>(f)[i];
                });
        case duckdb::LogicalTypeId::FLOAT:
            return readColumnInto(vec, count, targetRow, chunkOffset, dest,
                [](const duckdb::UnifiedVectorFormat& f, duckdb::idx_t i, double& out){
                    out = duckdb::UnifiedVectorFormat::GetData<float>(f)[i];
                });
        default:
            if(what) *what = "expected a floating-point column, got " + vec.GetType().ToString();
            return false;
    }
}

}

bool ReconstructedStructureContext::loadNeighborTopologyFromParquet(
    const std::string& neighborParquetPath,
    const LammpsParser::Frame& frame,
    StructureAnalysis& analysis,
    ReconstructedStructureContext& context,
    std::string* errorMessage
){
    if(neighborParquetPath.empty()){
        AnalysisDumpUtils::setError(errorMessage,
            "Neighbor topology Parquet path is empty (expected an upstream structure-identification "
            "stage to produce a 'neighbor_lattice' exposure).");
        return false;
    }

    const size_t natoms = static_cast<size_t>(frame.natoms);

    std::vector<std::vector<int>> indexStorage(MAX_NEIGHBORS, std::vector<int>(natoms, -1));
    std::array<std::array<std::vector<double>, MAX_NEIGHBORS>, 3> latticeStorage;
    for(int axis = 0; axis < 3; ++axis){
        for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
            latticeStorage[axis][slot].assign(natoms, 0.0);
        }
    }
    std::vector<std::int64_t> parquetIds(natoms, -1);

    try{
        auto db = Volt::Detail::openInMemoryDb();
        duckdb::Connection con(*db);

        const std::string sql =
            "SELECT * FROM read_parquet(" + sqlQuotePath(neighborParquetPath) + ")";
        auto result = con.Query(sql);
        if(result->HasError()){
            AnalysisDumpUtils::setError(errorMessage,
                "Failed to read neighbor topology Parquet: " + result->GetError());
            return false;
        }

        const auto rowCount = result->RowCount();
        if(rowCount != natoms){
            AnalysisDumpUtils::setError(errorMessage,
                "Neighbor topology Parquet row count (" + std::to_string(rowCount) +
                ") does not match frame atom count (" + std::to_string(natoms) + ").");
            return false;
        }

        const auto columnCount = result->ColumnCount();
        std::array<int, MAX_NEIGHBORS> indexCol;
        indexCol.fill(-1);
        std::array<std::array<int, MAX_NEIGHBORS>, 3> latticeCol;
        for(auto& axisCols : latticeCol) axisCols.fill(-1);
        int idCol = -1;

        int atomIndexCol = -1;
        const std::array<char, 3> axes = { 'x', 'y', 'z' };
        for(duckdb::idx_t c = 0; c < columnCount; ++c){
            const std::string& name = result->names[c];
            if(name == "id"){ idCol = static_cast<int>(c); continue; }
            if(name == "atom_index"){ atomIndexCol = static_cast<int>(c); continue; }
            for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
                if(name == AnalysisDumpUtils::neighborIndexName(slot)){
                    indexCol[slot] = static_cast<int>(c);
                }
                for(int axis = 0; axis < 3; ++axis){
                    if(name == AnalysisDumpUtils::neighborLatticeComponentName(axes[axis], slot)){
                        latticeCol[axis][slot] = static_cast<int>(c);
                    }
                }
            }
        }

        for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
            if(indexCol[slot] < 0){
                AnalysisDumpUtils::setError(errorMessage,
                    "Neighbor topology Parquet missing column '" +
                    AnalysisDumpUtils::neighborIndexName(slot) + "'");
                return false;
            }
            for(int axis = 0; axis < 3; ++axis){
                if(latticeCol[axis][slot] < 0){
                    AnalysisDumpUtils::setError(errorMessage,
                        "Neighbor topology Parquet missing column '" +
                        AnalysisDumpUtils::neighborLatticeComponentName(axes[axis], slot) + "'");
                    return false;
                }
            }
        }

        if(atomIndexCol < 0){
            AnalysisDumpUtils::setError(errorMessage,
                "Neighbor topology Parquet missing column 'atom_index'");
            return false;
        }

        std::vector<int> targetRow;
        std::string typeError;
        while(auto chunk = result->Fetch()){
            const duckdb::idx_t count = chunk->size();
            if(count == 0) continue;

            targetRow.assign(count, -1);
            {
                duckdb::UnifiedVectorFormat fmt;
                chunk->data[static_cast<duckdb::idx_t>(atomIndexCol)].ToUnifiedFormat(count, fmt);
                const auto typeId = chunk->data[static_cast<duckdb::idx_t>(atomIndexCol)].GetType().id();
                for(duckdb::idx_t r = 0; r < count; ++r){
                    const auto idx = fmt.sel->get_index(r);
                    if(!fmt.validity.RowIsValid(idx)) continue;
                    std::int64_t atomIndex = -1;
                    switch(typeId){
                        case duckdb::LogicalTypeId::TINYINT:
                            atomIndex = duckdb::UnifiedVectorFormat::GetData<int8_t>(fmt)[idx]; break;
                        case duckdb::LogicalTypeId::SMALLINT:
                            atomIndex = duckdb::UnifiedVectorFormat::GetData<int16_t>(fmt)[idx]; break;
                        case duckdb::LogicalTypeId::INTEGER:
                            atomIndex = duckdb::UnifiedVectorFormat::GetData<int32_t>(fmt)[idx]; break;
                        case duckdb::LogicalTypeId::BIGINT:
                            atomIndex = duckdb::UnifiedVectorFormat::GetData<int64_t>(fmt)[idx]; break;
                        case duckdb::LogicalTypeId::UTINYINT:
                            atomIndex = duckdb::UnifiedVectorFormat::GetData<uint8_t>(fmt)[idx]; break;
                        case duckdb::LogicalTypeId::USMALLINT:
                            atomIndex = duckdb::UnifiedVectorFormat::GetData<uint16_t>(fmt)[idx]; break;
                        case duckdb::LogicalTypeId::UINTEGER:
                            atomIndex = duckdb::UnifiedVectorFormat::GetData<uint32_t>(fmt)[idx]; break;
                        case duckdb::LogicalTypeId::UBIGINT:
                            atomIndex = static_cast<std::int64_t>(
                                duckdb::UnifiedVectorFormat::GetData<uint64_t>(fmt)[idx]); break;
                        default:
                            AnalysisDumpUtils::setError(errorMessage,
                                "Neighbor topology Parquet column 'atom_index' has unexpected type " +
                                chunk->data[static_cast<duckdb::idx_t>(atomIndexCol)].GetType().ToString());
                            return false;
                    }
                    if(atomIndex < 0 || static_cast<size_t>(atomIndex) >= natoms){
                        AnalysisDumpUtils::setError(errorMessage,
                            "Neighbor topology Parquet atom_index " + std::to_string(atomIndex) +
                            " is outside this frame's atom range (" + std::to_string(natoms) + ").");
                        return false;
                    }
                    targetRow[r] = static_cast<int>(atomIndex);
                }
            }

            if(idCol >= 0){
                if(!readIntegerColumn(chunk->data[static_cast<duckdb::idx_t>(idCol)], count,
                                    targetRow, 0, parquetIds.data(), &typeError)){
                    AnalysisDumpUtils::setError(errorMessage,
                        "Neighbor topology Parquet column 'id': " + typeError);
                    return false;
                }
            }
            for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
                if(!readIntegerColumn(chunk->data[static_cast<duckdb::idx_t>(indexCol[slot])], count,
                                  targetRow, 0, indexStorage[slot].data(), &typeError)){
                    AnalysisDumpUtils::setError(errorMessage,
                        "Neighbor topology Parquet column '" +
                        AnalysisDumpUtils::neighborIndexName(slot) + "': " + typeError);
                    return false;
                }
            }
            for(int axis = 0; axis < 3; ++axis){
                for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
                    if(!readDoubleColumn(chunk->data[static_cast<duckdb::idx_t>(latticeCol[axis][slot])],
                                         count, targetRow, 0,
                                         latticeStorage[axis][slot].data(), &typeError)){
                        AnalysisDumpUtils::setError(errorMessage,
                            "Neighbor topology Parquet column '" +
                            AnalysisDumpUtils::neighborLatticeComponentName(axes[axis], slot) +
                            "': " + typeError);
                        return false;
                    }
                }
            }
        }
    }catch(const std::exception& error){
        AnalysisDumpUtils::setError(errorMessage,
            std::string("Failed to read neighbor topology Parquet: ") + error.what());
        return false;
    }

    if(frame.ids.size() == natoms){
        for(size_t i = 0; i < natoms; ++i){
            if(parquetIds[i] >= 0 && parquetIds[i] != static_cast<std::int64_t>(frame.ids[i])){
                AnalysisDumpUtils::setError(errorMessage,
                    "Neighbor topology Parquet atom id mismatch at row " + std::to_string(i) +
                    " (parquet id " + std::to_string(parquetIds[i]) +
                    " != frame id " + std::to_string(frame.ids[i]) + ").");
                return false;
            }
        }
    }

    std::array<const int*, MAX_NEIGHBORS> neighborSlotPtrs{};
    for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
        neighborSlotPtrs[slot] = indexStorage[slot].data();
    }
    std::array<std::array<const double*, MAX_NEIGHBORS>, 3> latticePtrs{};
    for(int axis = 0; axis < 3; ++axis){
        for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
            latticePtrs[axis][slot] = latticeStorage[axis][slot].data();
        }
    }

    return assembleNeighborGraph(frame, analysis, context, neighborSlotPtrs, latticePtrs);
}

bool ReconstructedStructureLoader::load(
    const LammpsParser::Frame& frame,
    const std::string& neighborParquetPath,
    const ClusterGraphExportPaths& paths,
    StructureAnalysis& structureAnalysis,
    ReconstructedStructureContext& context,
    std::string* errorMessage
){
    auto mark = [](std::chrono::steady_clock::time_point& since){
        const auto now = std::chrono::steady_clock::now();
        const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - since).count();
        since = now;
        return ms;
    };
    auto since = std::chrono::steady_clock::now();

    if(!ReconstructedStructureContext::loadStructureAndClusterFromFrame(frame, context, errorMessage)){
        return false;
    }
    const auto fromFrameMs = mark(since);

    if(!ReconstructedStructureContext::loadNeighborTopologyFromParquet(
        neighborParquetPath, frame, structureAnalysis, context, errorMessage)){
        return false;
    }
    const auto neighborParquetMs = mark(since);

    if(!importClusterGraph(structureAnalysis, paths, errorMessage)){
        return false;
    }
    const auto clusterGraphMs = mark(since);

    spdlog::info(
        "  reconstruction: {}ms per-atom columns from dump, {}ms neighbor_lattice.parquet, {}ms cluster graph",
        fromFrameMs, neighborParquetMs, clusterGraphMs
    );

    return true;
}

}
