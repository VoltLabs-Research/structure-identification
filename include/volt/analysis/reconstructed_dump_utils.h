#pragma once

#include <memory>
#include <string>
#include <vector>

#include <volt/core/lammps_parser.h>
#include <volt/core/particle_property.h>
#include <volt/structures/crystal_structure_types.h>

namespace Volt::AnalysisDumpUtils{

std::string neighborIndexName(int slot);
std::string neighborLatticeComponentName(char axis, int slot);
std::vector<std::string> neighborIndexNames();
std::vector<std::string> neighborLatticeVectorNames();

void setError(std::string* errorMessage, const std::string& message);

bool tryParseInt(const std::string& text, int& value);
bool tryParseDouble(const std::string& text, double& value);

std::shared_ptr<ParticleProperty> makeIntProperty(std::size_t atomCount, int fillValue);
std::shared_ptr<ParticleProperty> makeIntProperty(const std::vector<int>& values);

bool extractIntegralColumn(
    const LammpsParser::Frame& frame,
    const std::string& columnName,
    std::vector<int>& output,
    std::string* errorMessage
);

bool extractDoubleColumn(
    const LammpsParser::Frame& frame,
    const std::string& columnName,
    std::vector<double>& output,
    std::string* errorMessage
);

}
