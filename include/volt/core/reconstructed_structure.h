#pragma once

#include <string>

#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cluster_graph_io.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/core/lammps_parser.h>

namespace Volt{

class ReconstructedStructureContext final : public StructureContext{
public:
    ReconstructedStructureContext(
        ParticleProperty* positions,
        const SimulationCell& cell
    );

    // Reads structure_type (optional) + cluster_id (required) from the dump frame.
    // Neighbor topology no longer rides in the dump — see loadNeighborTopologyFromParquet.
    static bool loadStructureAndClusterFromFrame(
        const LammpsParser::Frame& frame,
        ReconstructedStructureContext& context,
        std::string* errorMessage = nullptr
    );

    // Reads neighbor_indices_* + neighbor_lattice_* from the dedicated sidecar
    // Parquet (produced by streamNeighborTopologyToParquet) and rebuilds the CSR
    // neighbor graph + lattice-vector overrides. Atom positions still come from
    // `frame`; the parquet is joined by row order (ORDER BY atom_index) and
    // validated against frame.ids so the positional neighbor indices line up
    // with the consumer frame.
    static bool loadNeighborTopologyFromParquet(
        const std::string& neighborParquetPath,
        const LammpsParser::Frame& frame,
        StructureAnalysis& analysis,
        ReconstructedStructureContext& context,
        std::string* errorMessage = nullptr
    );

private:
    std::shared_ptr<ParticleProperty> _ownedStructureTypes;
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
