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

    static bool loadFromFrame(
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
        const ClusterGraphExportPaths& paths,
        StructureAnalysis& structureAnalysis,
        ReconstructedStructureContext& context,
        std::string* errorMessage = nullptr
    );
};

}
