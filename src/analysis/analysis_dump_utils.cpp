#include <volt/analysis/analysis_dump_utils.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

namespace Volt::AnalysisDumpUtils{

namespace{

template<typename IntType>
bool convertIntegral(double input, IntType& output){
    if(!std::isfinite(input)){
        return false;
    }

    const double rounded = std::round(input);
    if(std::abs(input - rounded) > 1e-6){
        return false;
    }

    const double minValue = static_cast<double>(std::numeric_limits<IntType>::lowest());
    const double maxValue = static_cast<double>(std::numeric_limits<IntType>::max());
    if(rounded < minValue || rounded > maxValue){
        return false;
    }

    output = static_cast<IntType>(rounded);
    return true;
}

}

std::string neighborIndexName(int slot){
    return "neighbor_indices_" + std::to_string(slot);
}

std::string neighborLatticeComponentName(char axis, int slot){
    return std::string("neighbor_lattice_") + axis + "_" + std::to_string(slot);
}

std::vector<std::string> neighborIndexNames(){
    std::vector<std::string> names;
    names.reserve(MAX_NEIGHBORS);
    for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
        names.push_back(neighborIndexName(slot));
    }
    return names;
}

std::vector<std::string> neighborLatticeVectorNames(){
    std::vector<std::string> names;
    names.reserve(MAX_NEIGHBORS * 3);
    for(int slot = 0; slot < MAX_NEIGHBORS; ++slot){
        names.push_back(neighborLatticeComponentName('x', slot));
        names.push_back(neighborLatticeComponentName('y', slot));
        names.push_back(neighborLatticeComponentName('z', slot));
    }
    return names;
}

void setError(std::string* errorMessage, const std::string& message){
    if(errorMessage){
        *errorMessage = message;
    }
}

bool tryParseInt(const std::string& text, int& value){
    try{
        value = std::stoi(text);
        return true;
    }catch(...){
        return false;
    }
}

bool tryParseDouble(const std::string& text, double& value){
    try{
        value = std::stod(text);
        return true;
    }catch(...){
        return false;
    }
}

std::shared_ptr<ParticleProperty> makeIntProperty(std::size_t atomCount, int fillValue){
    auto property = std::make_shared<ParticleProperty>(atomCount, DataType::Int, 1, 0, true);
    std::fill(property->dataInt(), property->dataInt() + atomCount, fillValue);
    return property;
}

std::shared_ptr<ParticleProperty> makeIntProperty(const std::vector<int>& values){
    auto property = std::make_shared<ParticleProperty>(values.size(), DataType::Int, 1, 0, true);
    std::copy(values.begin(), values.end(), property->dataInt());
    return property;
}

bool extractHeaderValue(
    const LammpsParser::Frame& frame,
    const std::string& headerName,
    std::string& value,
    std::string* errorMessage
){
    const std::string* header = frame.findHeaderProperty(headerName);
    if(!header){
        setError(errorMessage, "Missing required dump header '" + headerName + "'");
        return false;
    }

    value = *header;
    return true;
}

bool extractIntegralColumn(
    const LammpsParser::Frame& frame,
    const std::string& columnName,
    std::vector<int>& output,
    std::string* errorMessage
){
    const auto* column = frame.findAtomProperty(columnName);
    if(!column){
        setError(errorMessage, "Missing required dump column '" + columnName + "'");
        return false;
    }

    output.resize(static_cast<std::size_t>(frame.natoms));
    switch(column->dataType){
        case DataType::Int:
            if(column->ints.size() != output.size()){
                setError(errorMessage, "Column '" + columnName + "' size does not match atom count");
                return false;
            }
            std::copy(column->ints.begin(), column->ints.end(), output.begin());
            return true;
        case DataType::Int64:
            if(column->int64s.size() != output.size()){
                setError(errorMessage, "Column '" + columnName + "' size does not match atom count");
                return false;
            }
            for(std::size_t i = 0; i < output.size(); ++i){
                const std::int64_t value = column->int64s[i];
                if(value < static_cast<std::int64_t>(std::numeric_limits<int>::lowest()) ||
                   value > static_cast<std::int64_t>(std::numeric_limits<int>::max())){
                    setError(errorMessage, "Column '" + columnName + "' contains values out of range");
                    return false;
                }
                output[i] = static_cast<int>(value);
            }
            return true;
        case DataType::Double:
            if(column->doubles.size() != output.size()){
                setError(errorMessage, "Column '" + columnName + "' size does not match atom count");
                return false;
            }
            for(std::size_t i = 0; i < output.size(); ++i){
                if(!convertIntegral(column->doubles[i], output[i])){
                    setError(errorMessage, "Column '" + columnName + "' contains non-integral values");
                    return false;
                }
            }
            return true;
        case DataType::Void:
            break;
    }

    setError(errorMessage, "Unsupported data type for dump column '" + columnName + "'");
    return false;
}

bool extractDoubleColumn(
    const LammpsParser::Frame& frame,
    const std::string& columnName,
    std::vector<double>& output,
    std::string* errorMessage
){
    const auto* column = frame.findAtomProperty(columnName);
    if(!column){
        setError(errorMessage, "Missing required dump column '" + columnName + "'");
        return false;
    }

    output.resize(static_cast<std::size_t>(frame.natoms));
    switch(column->dataType){
        case DataType::Int:
            if(column->ints.size() != output.size()){
                setError(errorMessage, "Column '" + columnName + "' size does not match atom count");
                return false;
            }
            for(std::size_t i = 0; i < output.size(); ++i){
                output[i] = static_cast<double>(column->ints[i]);
            }
            return true;
        case DataType::Int64:
            if(column->int64s.size() != output.size()){
                setError(errorMessage, "Column '" + columnName + "' size does not match atom count");
                return false;
            }
            for(std::size_t i = 0; i < output.size(); ++i){
                output[i] = static_cast<double>(column->int64s[i]);
            }
            return true;
        case DataType::Double:
            if(column->doubles.size() != output.size()){
                setError(errorMessage, "Column '" + columnName + "' size does not match atom count");
                return false;
            }
            output = column->doubles;
            return true;
        case DataType::Void:
            break;
    }

    setError(errorMessage, "Unsupported data type for dump column '" + columnName + "'");
    return false;
}

}
