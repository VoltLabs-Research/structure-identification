#pragma once

#include <volt/analysis/analysis_context.h>
#include <volt/analysis/cluster_graph_export.h>
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
            nullptr,
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
    const AnalysisContext& context,
    StructureAnalysis& analysis,
    nlohmann::json& result,
    std::string* errorMessage
){
    if(outputBase.empty()){
        return true;
    }

    const std::string dumpPath = outputBase + "_annotated.dump";
    if(!context.writeDumpWithContext(frame, dumpPath, &analysis)){
        if(errorMessage){
            *errorMessage = "Failed to write " + dumpPath;
        }
        return false;
    }

    ClusterGraphExportPaths clusterGraphPaths;
    if(!exportClusterGraph(analysis.clusterGraph(), outputBase, &clusterGraphPaths)){
        if(errorMessage){
            *errorMessage = "Failed to export cluster graph tables";
        }
        return false;
    }

    result["annotated_dump"] = dumpPath;
    result["clusters_table"] = clusterGraphPaths.clustersTablePath;
    result["cluster_transitions_table"] = clusterGraphPaths.clusterTransitionsTablePath;
    return true;
}

}
