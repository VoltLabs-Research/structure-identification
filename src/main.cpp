#include <volt/cli/common.h>
#include <volt/core/frame_adapter.h>
#include <volt/analysis/analysis_context.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/utilities/json_utils.h>
#include <spdlog/spdlog.h>
#include <fstream>
#include <vector>
#include <memory>

using namespace Volt;
using namespace Volt::CLI;
using namespace Volt::Particles;

StructureAnalysis::Mode parseIdentificationMode(const std::string& value){
    if(value == "CNA") return StructureAnalysis::Mode::CNA;
    if(value == "PTM") return StructureAnalysis::Mode::PTM;
    if(value == "DIAMOND") return StructureAnalysis::Mode::DIAMOND;
    spdlog::warn("Unknown identification mode '{}', defaulting to CNA.", value);
    return StructureAnalysis::Mode::CNA;
}

void showUsage(const std::string& name){
    printUsageHeader(name, "Volt - Structure Identification");
    std::cerr
        << "  --mode <mode>     Identification mode. (CNA|PTM|DIAMOND) [default: CNA]\n"
        << "  --rmsd <float>    RMSD threshold for PTM. [default: 0.1]\n"
        << "  --threads <int>   Max worker threads (TBB/OMP). [default: auto]\n";
    printHelpOption();
}

bool exportPTMData(
    const AnalysisContext& context,
    const std::vector<int>& ids,
    const std::string& outputBase
){
    const auto ptmOrientation = context.ptmOrientation;
    const auto correspondences = context.correspondencesCode;

    if(!ptmOrientation || !correspondences){
        return true;
    }

    std::ofstream output(outputBase + "_ptm_data.msgpack", std::ios::binary);
    if(!output.is_open()){
        spdlog::error("Failed to open PTM output file.");
        return false;
    }

    MsgpackWriter writer(output);
    writer.write_map_header(1);
    writer.write_key("data");
    writer.write_array_header(JsonUtils::checked_u32_size(ids.size()));

    const bool includeStructureType = context.structureTypes && context.structureTypes->size() >= ids.size();
    for(size_t index = 0; index < ids.size(); ++index){
        writer.write_map_header(includeStructureType ? 4u : 3u);
        writer.write_key("id");
        writer.write_int(static_cast<int64_t>(ids[index]));

        writer.write_key("correspondence");
        writer.write_uint(static_cast<uint64_t>(correspondences->getInt64(index)));

        if(includeStructureType){
            writer.write_key("structure_type");
            writer.write_int(static_cast<int64_t>(context.structureTypes->getInt(static_cast<int>(index))));
        }

        writer.write_key("orientation");
        writer.write_array_header(4);
        for(int component = 0; component < 4; ++component){
            writer.write_double(ptmOrientation->getDoubleComponent(index, component));
        }
    }

    spdlog::info("PTM data written to {}_ptm_data.msgpack", outputBase);
    return true;
}

int main(int argc, char* argv[]){
    if(argc < 2){
        showUsage(argv[0]);
        return 1;
    }

    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);

    if(hasOption(opts, "--help") || filename.empty()){
        showUsage(argv[0]);
        return filename.empty() ? 1 : 0;
    }

    auto parallel = initParallelism(opts, false);
    initLogging("volt-structure-identification", parallel.threads);

    LammpsParser::Frame frame;
    if(!parseFrame(filename, frame)) return 1;

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions){
        spdlog::error("Failed to create position property.");
        return 1;
    }

    if(!FrameAdapter::validateSimulationCell(frame.simulationCell)){
        spdlog::error("Invalid simulation cell.");
        return 1;
    }

    auto structureTypes = std::make_unique<ParticleProperty>(frame.natoms, DataType::Int, 1, 0, true);
    std::vector<Matrix3> preferredOrientations;
    preferredOrientations.push_back(Matrix3::Identity());

    AnalysisContext context(
        positions.get(),
        frame.simulationCell,
        LATTICE_FCC,
        nullptr,
        structureTypes.get(),
        std::move(preferredOrientations)
    );

    const auto mode = parseIdentificationMode(getString(opts, "--mode", "CNA"));
    const float rmsd = static_cast<float>(getDouble(opts, "--rmsd", 0.1f));

    spdlog::info("Starting structure identification...");
    StructureAnalysis analysis(context, true, mode, rmsd);
    analysis.identifyStructures();

    std::vector<int> atomStructureTypes(static_cast<size_t>(frame.natoms), static_cast<int>(StructureType::OTHER));
    for(int atomIndex = 0; atomIndex < frame.natoms; ++atomIndex){
        atomStructureTypes[static_cast<size_t>(atomIndex)] = context.structureTypes->getInt(atomIndex);
    }

    if(!JsonUtils::writeJsonMsgpackToFile(
        analysis.getAtomsData(frame, &atomStructureTypes),
        outputBase + "_atoms.msgpack",
        false
    )){
        spdlog::error("Failed to write atoms.msgpack");
        return 1;
    }

    if(!JsonUtils::writeJsonMsgpackToFile(
        analysis.getStructureStatisticsJson(),
        outputBase + "_structure_analysis_stats.msgpack",
        false
    )){
        spdlog::error("Failed to write structure_analysis_stats.msgpack");
        return 1;
    }

    std::vector<int> atomIds(static_cast<size_t>(frame.natoms));
    for(int atomIndex = 0; atomIndex < frame.natoms; ++atomIndex){
        atomIds[static_cast<size_t>(atomIndex)] = atomIndex < static_cast<int>(frame.ids.size())
            ? frame.ids[static_cast<size_t>(atomIndex)]
            : atomIndex;
    }

    if(mode == StructureAnalysis::Mode::PTM){
        if(!exportPTMData(context, atomIds, outputBase)){
            return 1;
        }
    }

    spdlog::info("Structure identification completed.");
    return 0;
}
