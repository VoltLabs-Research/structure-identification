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

json buildPtmPerAtomProperties(
    const AnalysisContext& context,
    const LammpsParser::Frame& frame,
    const std::vector<int>& structureTypes
){
    const auto ptmOrientation = context.ptmOrientation;
    const auto correspondences = context.correspondencesCode;

    json perAtom = json::array();

    for(size_t index = 0; index < frame.natoms; ++index){
        int atomId = index < frame.ids.size()
            ? frame.ids[index]
            : static_cast<int>(index);

        json atom;
        atom["id"] = atomId;
        atom["structure_type"] = structureTypes[index];

        if(index < static_cast<size_t>(frame.positions.size())){
            const auto& pos = frame.positions[index];
            atom["pos"] = {pos.x(), pos.y(), pos.z()};
        }else{
            atom["pos"] = {0.0, 0.0, 0.0};
        }

        if(correspondences){
            atom["correspondence"] = static_cast<uint64_t>(correspondences->getInt64(index));
        }

        if(ptmOrientation){
            atom["orientation"] = {
                ptmOrientation->getDoubleComponent(index, 0),
                ptmOrientation->getDoubleComponent(index, 1),
                ptmOrientation->getDoubleComponent(index, 2),
                ptmOrientation->getDoubleComponent(index, 3)
            };
        }

        perAtom.push_back(atom);
    }

    return perAtom;
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

    json result;
    result["main_listing"] = analysis.buildMainListing();

    if(mode == StructureAnalysis::Mode::PTM){
        result["per-atom-properties"] = buildPtmPerAtomProperties(context, frame, atomStructureTypes);
    }else{
        result["per-atom-properties"] = analysis.getPerAtomProperties(frame, &atomStructureTypes);
    }

    const std::string outputPath = outputBase + "_structure_identification.msgpack";
    if(!JsonUtils::writeJsonMsgpackToFile(result, outputPath, false)){
        spdlog::error("Failed to write {}", outputPath);
        return 1;
    }

    spdlog::info("Structure identification completed. Output: {}", outputPath);
    return 0;
}
