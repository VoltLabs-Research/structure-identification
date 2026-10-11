#pragma once
#include <volt/core/particle_property.h>
#include <volt/core/simulation_cell.h>
#include <volt/math/matrix3.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/core/lammps_parser.h>

#include <memory>
#include <string>
#include <vector>

namespace Volt{

class StructureAnalysis;

class StructureContext{
public:
    ParticleProperty* positions;
    ParticleProperty* structureTypes;

    double maximumNeighborDistance;
    std::shared_ptr<ParticleProperty> neighborOffsets;
    std::shared_ptr<ParticleProperty> neighborIndices;
    std::shared_ptr<ParticleProperty> neighborCounts;
    std::shared_ptr<ParticleProperty> atomClusters;
    std::shared_ptr<ParticleProperty> atomAllowedSymmetryMasks;
    std::shared_ptr<ParticleProperty> atomSymmetryPermutations;

    const SimulationCell& simCell;
    LatticeStructureType inputCrystalType;
    std::vector<Matrix3> preferredCrystalOrientations;

    StructureContext(
        ParticleProperty* pos,
        const SimulationCell& cell,
        LatticeStructureType crystalType = LATTICE_OTHER,
        ParticleProperty* outputStructures = nullptr,
        std::vector<Matrix3>&& preferredOrientations = {}
    );

    virtual ~StructureContext() = default;

    size_t atomCount() const{
        return positions->size();
    }
};

class AnalysisContext : public StructureContext{
public:
    struct ExportedContext{
        std::vector<LammpsParser::ExtraHeader> headers;
        std::vector<LammpsParser::ExtraColumn> columns;
        std::vector<std::shared_ptr<ParticleProperty>> ownedProperties;
    };

    struct ExtraScalarColumn{
        std::string name;
        std::shared_ptr<ParticleProperty> property;
    };

    ExportedContext exportContext(
        const StructureAnalysis* analysis = nullptr,
        const std::vector<ExtraScalarColumn>& extraColumns = {}
    ) const;
    bool writeDumpWithContext(
        const LammpsParser::Frame& frame,
        const std::string& outputFilename,
        const StructureAnalysis* analysis = nullptr,
        const std::vector<ExtraScalarColumn>& extraColumns = {},
        const std::vector<bool>* atomMask = nullptr
    ) const;

    AnalysisContext(
        ParticleProperty* pos,
        const SimulationCell& cell,
        LatticeStructureType crystalType,
        ParticleProperty* outputStructures,
        std::vector<Matrix3>&& preferredOrientations
    );
};

int contractAtomId(const LammpsParser::Frame& frame, std::size_t atomIndex);

LammpsParser::Frame maskFrame(const LammpsParser::Frame& frame, const std::vector<bool>& atomMask);

bool streamNeighborTopologyToParquet(
    const std::string& filePath,
    const LammpsParser::Frame& frame,
    const AnalysisContext& context,
    const StructureAnalysis& analysis,
    const std::vector<bool>* atomMask = nullptr
);

}
