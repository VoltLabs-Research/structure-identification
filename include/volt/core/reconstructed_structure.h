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
