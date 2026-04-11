#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/reconstructed_dump_utils.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/structures/crystal_structure_types.h>

#include <algorithm>
#include <stdexcept>
#include <string>

namespace Volt{

namespace AnalysisContextDetail{

std::shared_ptr<ParticleProperty> makeNeighborIndicesProperty(const AnalysisContext& context){
    const std::size_t atomCount = context.atomCount();
    auto property = std::make_shared<ParticleProperty>(
        atomCount,
        DataType::Int,
        MAX_NEIGHBORS,
        0,
        true
    );
    std::fill(
        property->dataInt(),
        property->dataInt() + atomCount * MAX_NEIGHBORS,
        -1
    );

    if(!context.neighborCounts || !context.neighborOffsets || !context.neighborIndices){
        return property;
    }

    const int* counts = context.neighborCounts->constDataInt();
    const int* offsets = context.neighborOffsets->constDataInt();
    const int* indices = context.neighborIndices->constDataInt();
    const std::size_t totalNeighborEntries = context.neighborIndices->size();
    int* destination = property->dataInt();

    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        const int count = std::max(0, std::min(counts[atomIndex], static_cast<int>(MAX_NEIGHBORS)));
        const int start = offsets[atomIndex];
        if(start < 0){
            continue;
        }
        for(int neighborSlot = 0; neighborSlot < count; ++neighborSlot){
            const std::size_t sourceIndex = static_cast<std::size_t>(start + neighborSlot);
            if(sourceIndex >= totalNeighborEntries){
                break;
            }
            destination[atomIndex * MAX_NEIGHBORS + neighborSlot] = indices[sourceIndex];
        }
    }

    return property;
}

std::shared_ptr<ParticleProperty> makeNeighborLatticeVectorProperty(
    const AnalysisContext& context,
    const StructureAnalysis& analysis
){
    auto property = std::make_shared<ParticleProperty>(
        context.atomCount(),
        DataType::Double,
        MAX_NEIGHBORS * 3,
        0.0,
        true
    );

    for(std::size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
        const int structureType = context.structureTypes->getInt(atomIndex);
        if(structureType == LATTICE_OTHER){
            continue;
        }

        const int neighborCount = context.neighborCounts
            ? context.neighborCounts->getInt(static_cast<int>(atomIndex))
            : 0;
        const int maxCount = analysis.hasNeighborLatticeVectorOverrides()
            ? static_cast<int>(MAX_NEIGHBORS)
            : analysis.coordinationNumber(structureType);
        const int exportableCount = std::max(0, std::min(neighborCount, maxCount));

        for(int neighborSlot = 0; neighborSlot < exportableCount; ++neighborSlot){
            const Vector3& vector = analysis.neighborLatticeVector(static_cast<int>(atomIndex), neighborSlot);
            const std::size_t baseComponent = static_cast<std::size_t>(neighborSlot) * 3;
            property->setDoubleComponent(atomIndex, baseComponent + 0, vector.x());
            property->setDoubleComponent(atomIndex, baseComponent + 1, vector.y());
            property->setDoubleComponent(atomIndex, baseComponent + 2, vector.z());
        }
    }

    return property;
}

}

StructureContext::StructureContext(
    ParticleProperty* pos,
    const SimulationCell& cell,
    LatticeStructureType crystalType,
    ParticleProperty* selection,
    ParticleProperty* outputStructures,
    std::vector<Matrix3>&& preferredOrientations
) : 
    positions(pos),
    structureTypes(outputStructures),
    particleSelection(selection),
    maximumNeighborDistance(0.0),
    simCell(cell),
    inputCrystalType(crystalType),
    preferredCrystalOrientations(std::move(preferredOrientations))
{
    if(!positions){
        throw std::invalid_argument("Invalid positions");
    }
}

using namespace AnalysisContextDetail;

AnalysisContext::AnalysisContext(
    ParticleProperty* pos,
    const SimulationCell& cell,
    LatticeStructureType crystalType,
    ParticleProperty* selection,
    ParticleProperty* outputStructures,
    std::vector<Matrix3>&& preferredOrientations
) :
    StructureContext(
        pos,
        cell,
        crystalType,
        selection,
        outputStructures,
        std::move(preferredOrientations)
    )
{
    if(!structureTypes){
        throw std::invalid_argument("Invalid structure types");
    }
    const size_t numAtoms = atomCount();
    atomClusters = std::make_shared<ParticleProperty>(numAtoms, DataType::Int, 1, 0, true);
    atomAllowedSymmetryMasks = std::make_shared<ParticleProperty>(numAtoms, DataType::Int64, 1, 0, true);
    atomSymmetryPermutations = std::make_shared<ParticleProperty>(numAtoms, DataType::Int, 1, 0, true);
    std::fill(
        atomSymmetryPermutations->dataInt(),
        atomSymmetryPermutations->dataInt() + numAtoms,
        -1
    );

    if(numAtoms > 0){
        std::fill(
            structureTypes->dataInt(),
            structureTypes->dataInt() + numAtoms,
            LATTICE_OTHER
        );
    }
}

AnalysisContext::ExportedContext AnalysisContext::exportContext(const StructureAnalysis* analysis) const{
    ExportedContext exported;
    auto& columns = exported.columns;
    auto& ownedProperties = exported.ownedProperties;
    const std::size_t count = atomCount();

    auto writeIntColumn = [&](const char* name, const std::shared_ptr<ParticleProperty>& property, int fillValue){
        auto exportedProperty = property ? property : AnalysisDumpUtils::makeIntProperty(count, fillValue);
        if(!property){
            ownedProperties.push_back(exportedProperty);
        }
        LammpsParser::writeColumn(columns, { name }, exportedProperty);
    };
    auto writeRawIntColumn = [&](const char* name, const ParticleProperty* property, int fillValue){
        if(!property){
            writeIntColumn(name, nullptr, fillValue);
            return;
        }

        auto exportedProperty = std::make_shared<ParticleProperty>(count, DataType::Int, 1, 0, true);
        std::copy(
            property->constDataInt(),
            property->constDataInt() + count,
            exportedProperty->dataInt()
        );
        ownedProperties.push_back(exportedProperty);
        LammpsParser::writeColumn(columns, { name }, exportedProperty);
    };
    writeRawIntColumn("structure_type", structureTypes, LATTICE_OTHER);
    writeIntColumn("cluster_id", atomClusters, 0);

    auto neighborIndicesProperty = makeNeighborIndicesProperty(*this);
    ownedProperties.push_back(neighborIndicesProperty);
    LammpsParser::writeColumn(columns, AnalysisDumpUtils::neighborIndexNames(), neighborIndicesProperty);

    if(analysis){
        auto neighborLatticeVectors = makeNeighborLatticeVectorProperty(*this, *analysis);
        ownedProperties.push_back(neighborLatticeVectors);
        LammpsParser::writeColumn(columns, AnalysisDumpUtils::neighborLatticeVectorNames(), neighborLatticeVectors);
    }

    return exported;
}

bool AnalysisContext::writeDumpWithContext(
    const LammpsParser::Frame& frame,
    const std::string& outputFilename,
    const StructureAnalysis* analysis
) const{
    LammpsParser parser;
    const auto exported = exportContext(analysis);

    const std::vector<int>* atomIds = nullptr;
    std::vector<int> generatedAtomIds;

    if(frame.ids.size() == atomCount()){
        atomIds = &frame.ids;
    }else{
        generatedAtomIds.resize(atomCount());
        for(std::size_t atomIndex = 0; atomIndex < atomCount(); ++atomIndex){
            generatedAtomIds[atomIndex] = static_cast<int>(atomIndex);
        }
        atomIds = &generatedAtomIds;
    }

    return parser.writeFileMergedWithExtraColumns(
        outputFilename,
        frame,
        *atomIds,
        exported.columns,
        exported.headers,
        true
    );
}

}
