#include <volt/core/reconstructed_structure.h>
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

namespace Volt{

namespace{

// Builds the CSR neighbor graph + lattice-vector overrides from raw per-slot
// column pointers. `neighborSlotPtrs[s][atomIndex]` is the positional index of
// atom `atomIndex`'s s-th neighbor (negative ⇒ no neighbor in that slot);
// `latticePtrs[axis][s][atomIndex]` is that neighbor's ideal lattice vector
// component. Positions for the maximum-neighbor-distance scan come from `frame`.
// This is the exact assembly that previously lived inline in loadFromFrame —
// the only change is its inputs now come from the sidecar Parquet, not the dump.
bool assembleNeighborGraph(
    const LammpsParser::Frame& frame,
    StructureAnalysis& analysis,
    ReconstructedStructureContext& context,
    const std::array<const int*, MAX_NEIGHBORS>& neighborSlotPtrs,
    const std::array<std::array<const double*, MAX_NEIGHBORS>, 3>& latticePtrs
){
    const size_t natoms = static_cast<size_t>(frame.natoms);

    // Parallel neighbor counting
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

    // Prefix sum for offsets
    std::vector<int> compactOffsets(natoms + 1, 0);
    int totalNeighborEntries = 0;
    for(size_t i = 0; i < natoms; ++i){
        compactOffsets[i] = totalNeighborEntries;
        totalNeighborEntries += neighborCounts[i];
    }
    compactOffsets[natoms] = totalNeighborEntries;

    context.neighborCounts = AnalysisDumpUtils::makeIntProperty(neighborCounts);
    context.neighborOffsets = AnalysisDumpUtils::makeIntProperty(compactOffsets);

    // Parallel compact neighbor indices
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

    // Parallel max neighbor distance
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

    // Parallel vector overrides assembly
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

    // Materialise the 18 index columns + 54 lattice columns into stable per-column
    // buffers, ordered by atom_index so row i corresponds to frame atom i.
    std::vector<std::vector<int>> indexStorage(MAX_NEIGHBORS, std::vector<int>(natoms, -1));
    std::array<std::array<std::vector<double>, MAX_NEIGHBORS>, 3> latticeStorage;
    for(int axis = 0; axis < 3; ++axis){
        for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
            latticeStorage[axis][slot].assign(natoms, 0.0);
        }
    }
    std::vector<std::int64_t> parquetIds(natoms, -1);

    try{
        duckdb::DuckDB db(nullptr);
        duckdb::Connection con(db);

        const std::string sql =
            "SELECT * FROM read_parquet(" + sqlQuotePath(neighborParquetPath) + ") ORDER BY atom_index";
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

        // Map column name → result column index.
        const auto columnCount = result->ColumnCount();
        std::array<int, MAX_NEIGHBORS> indexCol;
        indexCol.fill(-1);
        std::array<std::array<int, MAX_NEIGHBORS>, 3> latticeCol;
        for(auto& axisCols : latticeCol) axisCols.fill(-1);
        int idCol = -1;

        const std::array<char, 3> axes = { 'x', 'y', 'z' };
        for(duckdb::idx_t c = 0; c < columnCount; ++c){
            const std::string& name = result->names[c];
            if(name == "id"){ idCol = static_cast<int>(c); continue; }
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

        for(duckdb::idx_t row = 0; row < rowCount; ++row){
            if(idCol >= 0){
                const duckdb::Value v = result->GetValue(static_cast<duckdb::idx_t>(idCol), row);
                parquetIds[row] = v.IsNull() ? -1 : v.GetValue<std::int64_t>();
            }
            for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
                const duckdb::Value v = result->GetValue(static_cast<duckdb::idx_t>(indexCol[slot]), row);
                indexStorage[slot][row] = v.IsNull() ? -1 : static_cast<int>(v.GetValue<std::int64_t>());
            }
            for(int axis = 0; axis < 3; ++axis){
                for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
                    const duckdb::Value v = result->GetValue(static_cast<duckdb::idx_t>(latticeCol[axis][slot]), row);
                    latticeStorage[axis][slot][row] = v.IsNull() ? 0.0 : v.GetValue<double>();
                }
            }
        }
    }catch(const std::exception& error){
        AnalysisDumpUtils::setError(errorMessage,
            std::string("Failed to read neighbor topology Parquet: ") + error.what());
        return false;
    }

    // Safety: positional neighbor indices are only valid if the Parquet rows line
    // up with this frame's atoms. atom_index ordering guarantees row i == atom i;
    // additionally assert the id columns match so a stale/mismatched sidecar fails
    // loudly instead of producing wrong dislocations.
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
    if(!ReconstructedStructureContext::loadStructureAndClusterFromFrame(frame, context, errorMessage)){
        return false;
    }

    if(!ReconstructedStructureContext::loadNeighborTopologyFromParquet(
        neighborParquetPath, frame, structureAnalysis, context, errorMessage)){
        return false;
    }

    if(!importClusterGraph(structureAnalysis, paths, errorMessage)){
        return false;
    }

    return true;
}

}
