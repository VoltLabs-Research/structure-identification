#pragma once

#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cluster_graph_io.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/core/frame_adapter.h>

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

inline bool appendClusterOutputs(
    const LammpsParser::Frame& frame,
    const std::string& outputBase,
    const std::string& inputDumpPath,
    const AnalysisContext& context,
    StructureAnalysis& analysis,
    nlohmann::json& result,
    std::string* errorMessage,
    const std::vector<AnalysisContext::ExtraScalarColumn>& extraColumns = {}
){
    if(!context.writeDumpWithContext(frame, inputDumpPath, &analysis, extraColumns)){
        if(errorMessage){
            *errorMessage = "Failed to write " + inputDumpPath;
        }
        return false;
    }
    result["input_dump"] = inputDumpPath;

    if(outputBase.empty()){
        return true;
    }

    ClusterGraphExportPaths clusterGraphPaths;
    if(!exportClusterGraph(analysis.clusterGraph(), outputBase, &clusterGraphPaths)){
        if(errorMessage){
            *errorMessage = "Failed to export cluster graph tables";
        }
        return false;
    }

    result["clusters_table"] = clusterGraphPaths.clustersTablePath;
    result["cluster_transitions_table"] = clusterGraphPaths.clusterTransitionsTablePath;

    // Neighbor topology sidecar: the per-atom neighbor graph + ideal lattice
    // vectors that OpenDXA / ElasticStrain / LineReconstructionDXA consume via
    // the `neighbor_lattice` inferFromContext exposure (kept out of atoms.parquet
    // and the annotated dump).
    const std::string neighborLatticePath = outputBase + "_neighbor_lattice.parquet";
    if(!streamNeighborTopologyToParquet(neighborLatticePath, frame, context, analysis)){
        if(errorMessage){
            *errorMessage = "Failed to write " + neighborLatticePath;
        }
        return false;
    }
    result["neighbor_lattice"] = neighborLatticePath;
    return true;
}

}
