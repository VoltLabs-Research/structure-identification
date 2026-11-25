#pragma once

#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cluster_graph_io.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/progress.h>

#include <memory>
#include <string>
#include <utility>

#include <nlohmann/json.hpp>

namespace Volt::AnalysisPipelineUtils{

struct PreparedAnalysisSession{
    std::shared_ptr<ParticleProperty> positions;
    std::unique_ptr<ParticleProperty> structureTypes;
    AnalysisContext context;

    PreparedAnalysisSession(
        const LammpsParser::Frame& frame,
        std::shared_ptr<ParticleProperty> preparedPositions,
        LatticeStructureType inputCrystalStructure
    )
        : positions(std::move(preparedPositions))
        , structureTypes(std::make_unique<ParticleProperty>(frame.natoms, DataType::Int, 1, 0, true))
        , context(
            positions.get(),
            frame.simulationCell,
            inputCrystalStructure,
            structureTypes.get(),
            std::vector<Matrix3>{Matrix3::Identity()}
        ){}
};

inline std::unique_ptr<PreparedAnalysisSession> prepareAnalysisSession(
    const LammpsParser::Frame& frame,
    LatticeStructureType inputCrystalStructure,
    std::string* errorMessage
){
    FrameAdapter::PreparedAnalysisInput prepared;
    if(!FrameAdapter::prepareAnalysisInput(frame, prepared, errorMessage)){
        return nullptr;
    }

    return std::make_unique<PreparedAnalysisSession>(
        frame,
        std::move(prepared.positions),
        inputCrystalStructure
    );
}

// The files DXA reads an analysis back from, as appendClusterOutputs() wrote them: the annotated dump and, next to
// an output base, the cluster graph tables and the neighbor lattice. False when one of them could not be written.
struct ClusterOutputs{
    std::string annotatedDump;
    ClusterGraphExportPaths clusterGraph;
    std::string neighborLattice;
    bool written = false;

    explicit operator bool() const{ return written; }
};

inline ClusterOutputs appendClusterOutputs(
    const LammpsParser::Frame& frame,
    const std::string& outputBase,
    const std::string& inputDumpPath,
    const AnalysisContext& context,
    StructureAnalysis& analysis,
    nlohmann::json& result,
    std::string* errorMessage,
    const std::vector<AnalysisContext::ExtraScalarColumn>& extraColumns = {},
    const std::vector<bool>* atomMask = nullptr
){
    ProgressSteps steps({
        {"Writing the annotated dump", 1.0},
        {"Writing the cluster graph tables", 0.2},
        {"Writing the neighbor lattice", 4.0},
    });
    ClusterOutputs outputs;
    if(!context.writeDumpWithContext(frame, inputDumpPath, &analysis, extraColumns, atomMask)){
        if(errorMessage){
            *errorMessage = "Failed to write " + inputDumpPath;
        }
        return outputs;
    }
    outputs.annotatedDump = inputDumpPath;
    result["input_dump"] = inputDumpPath;

    if(outputBase.empty()){
        outputs.written = true;
        return outputs;
    }
    steps.next();

    if(!exportClusterGraph(analysis.clusterGraph(), outputBase, &outputs.clusterGraph)){
        if(errorMessage){
            *errorMessage = "Failed to export cluster graph tables";
        }
        return outputs;
    }

    result["clusters_table"] = outputs.clusterGraph.clustersTablePath;
    result["cluster_transitions_table"] = outputs.clusterGraph.clusterTransitionsTablePath;
    steps.next();

    const std::string neighborLatticePath = outputBase + "_neighbor_lattice.parquet";
    if(!streamNeighborTopologyToParquet(neighborLatticePath, frame, context, analysis, atomMask, errorMessage)){
        return outputs;
    }
    outputs.neighborLattice = neighborLatticePath;
    result["neighbor_lattice"] = neighborLatticePath;
    outputs.written = true;
    return outputs;
}

}
