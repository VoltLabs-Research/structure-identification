#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cluster_graph_io.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/core/lammps_parser.h>

namespace duckdb{
class Vector;
}

namespace Volt{

class ReconstructedStructureContext final : public StructureContext{
public:
    ReconstructedStructureContext(
        ParticleProperty* positions,
        const SimulationCell& cell
    );

    static bool loadStructureAndClusterFromFrame(
        const LammpsParser::Frame& frame,
        ReconstructedStructureContext& context,
        std::string* errorMessage = nullptr
    );

    static bool loadNeighborTopologyFromParquet(
        const std::string& neighborParquetPath,
        const LammpsParser::Frame& frame,
        StructureAnalysis& analysis,
        ReconstructedStructureContext& context,
        std::string* errorMessage = nullptr
    );

private:
    static bool assembleNeighborGraph(
        const LammpsParser::Frame& frame,
        StructureAnalysis& analysis,
        ReconstructedStructureContext& context,
        const std::array<const int*, MAX_NEIGHBORS>& neighborSlotPtrs,
        const std::array<std::array<const double*, MAX_NEIGHBORS>, 3>& latticePtrs
    );

    // Each read stores the valid rows of a parquet column chunk at the atoms targetRow gives them (none for -1).
    template<typename Dest, typename Convert>
    static bool readColumnInto(duckdb::Vector& vec, std::uint64_t count, const std::vector<int>& targetRow,
                               std::uint64_t chunkOffset, Dest* dest, Convert convert);
    template<typename Dest>
    static bool readIntegerColumn(duckdb::Vector& vec, std::uint64_t count, const std::vector<int>& targetRow,
                                  std::uint64_t chunkOffset, Dest* dest, std::string* what);
    static bool readDoubleColumn(duckdb::Vector& vec, std::uint64_t count, const std::vector<int>& targetRow,
                                 std::uint64_t chunkOffset, double* dest, std::string* what);

    std::shared_ptr<ParticleProperty> _ownedStructureTypes;
};

// A cluster graph that clusters and transitions are imported into, with the ones it held before the import.
struct ClusterGraphImport{
    const ClusterGraph& graph;
    std::size_t clustersBefore = 0;
    std::size_t transitionsBefore = 0;
};

class ReconstructedStructureLoader{
public:
    static bool load(
        const LammpsParser::Frame& frame,
        const std::string& neighborParquetPath,
        const ClusterGraphExportPaths& paths,
        StructureAnalysis& structureAnalysis,
        ReconstructedStructureContext& context,
        std::string* errorMessage = nullptr
    );
};

}
