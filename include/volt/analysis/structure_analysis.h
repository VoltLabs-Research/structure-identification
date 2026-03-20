#pragma once

#include <volt/core/volt.h>
#include <volt/core/particle_property.h>
#include <volt/core/simulation_cell.h>
#include <volt/structures/cluster_graph.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/analysis/analysis_context.h>
#include <volt/core/lammps_parser.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace Volt{

class StructureAnalysis{
public:
	enum Mode{
		CNA,
		PTM,
		DIAMOND,
	};

	StructureAnalysis(
		AnalysisContext& context,
		bool identifyPlanarDefects, 
		Mode identificationMode,
		float rmsd
	);
	~StructureAnalysis();

	void identifyStructures();

	json getPerAtomProperties(
		const LammpsParser::Frame &frame,
		const std::vector<int>* structureTypes
	);

	void identifyStructuresCNA();
	void computeMaximumNeighborDistanceFromPTM();
	void determineLocalStructuresWithPTM();
	void appendNeighbors(const std::vector<std::vector<int>>& extras);

	int numberOfNeighbors(int atomIndex) const {
		assert(_context.neighborCounts);
		return _context.neighborCounts->getInt(atomIndex);
	}
	
	int getNeighbor(int centralAtomIndex, int neighborListIndex) const{
		assert(_context.neighborOffsets && _context.neighborIndices);
		const int count = _context.neighborCounts->getInt(centralAtomIndex);
		if(neighborListIndex < 0 || neighborListIndex >= count){
			return -1;
		}
		const int* offsets = _context.neighborOffsets->constDataInt();
		const int* indices = _context.neighborIndices->constDataInt();
		const int start = offsets[centralAtomIndex];
		return indices[start + neighborListIndex];
	}

	int findNeighbor(int centralAtomIndex, int neighborAtomIndex) const{
		assert(_context.neighborOffsets && _context.neighborIndices);
		const int count = _context.neighborCounts->getInt(centralAtomIndex);
		const int* offsets = _context.neighborOffsets->constDataInt();
		const int* indices = _context.neighborIndices->constDataInt();
		const int start = offsets[centralAtomIndex];
		for(int index = 0; index < count; index++){
			if(indices[start + index] == neighborAtomIndex){
				return index;
			}
		}
		return -1;
	}

	double maximumNeighborDistance() const{
		return _maximumNeighborDistance;
	}

	bool usingPTM() const{
		return _identificationMode == StructureAnalysis::Mode::PTM;
	}

	const AnalysisContext& context() const{
		return _context;
	}

	const ClusterGraph& clusterGraph() const{
		return *_clusterGraph;
	}

	ClusterGraph& clusterGraph(){
		return *_clusterGraph;
	}

	Cluster* atomCluster(int atomIndex) const{
		return clusterGraph().findCluster(_context.atomClusters->getInt(atomIndex));
	}
	
	void freeNeighborLists(){
		_context.neighborOffsets.reset();
		_context.neighborIndices.reset();
		_context.neighborCounts.reset();
		_context.atomSymmetryPermutations.reset();
	}

	void freePTMData(){
		_context.ptmRmsd.reset();
		_context.ptmOrientation.reset();
		_context.ptmDeformationGradient.reset();
		_context.correspondencesCode.reset();
		_context.templateIndex.reset();
	}

	void setIdentificationMode(Mode identificationMode){
		_identificationMode = identificationMode;
	}

	int findClosestSymmetryPermutation(int structureType, const Matrix3& rotation);
	int coordinationNumber(int structureType) const;
	int commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const;
	int symmetryPermutationCount(int structureType) const;
	int symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const;
	const Matrix3& symmetryTransformation(int structureType, int symmetryIndex) const;
	int symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const;
	const Vector3& latticeVector(int structureType, int latticeVectorIndex) const;

	// Returns the ideal lattice vector associated with a neighbor bond
	const Vector3& neighborLatticeVector(int centralAtomIndex, int neighborIndex) const;

	void calculateStructureStatistics() const {
        _structureStatistics.clear();
        
        const size_t N = _context.atomCount();
        for (size_t i = 0; i < N; ++i) {
            int structureType = _context.structureTypes->getInt(i);
            _structureStatistics[structureType]++;
        }
        
        _statisticsValid = true;
    }
    
    const std::map<int, int>& getStructureStatistics() const {
        if (!_statisticsValid) {
            calculateStructureStatistics();
        }
        return _structureStatistics;
    }
    
    std::map<std::string, int> getNamedStructureStatistics() const {
        if (!_statisticsValid) {
            calculateStructureStatistics();
        }
        
        std::map<std::string, int> namedStats;
        
        for (const auto& [structureType, count] : _structureStatistics) {
            std::string name = getStructureTypeName(structureType);
            namedStats[name] = count;
        }
        
        return namedStats;
    }
    
	void invalidateStatistics() {
		_statisticsValid = false;
	}
    
	std::string getStructureTypeName(int structureType) const {
		if (usingPTM()) {
            switch (static_cast<StructureType>(structureType)) {
                case StructureType::OTHER: return "OTHER";
                case StructureType::FCC: return "FCC";
                case StructureType::HCP: return "HCP";
                case StructureType::BCC: return "BCC";
                case StructureType::ICO: return "ICO";
                case StructureType::SC: return "SC";
                case StructureType::CUBIC_DIAMOND: return "CUBIC_DIAMOND";
				case StructureType::CUBIC_DIAMOND_FIRST_NEIGH: return "CUBIC_DIAMOND_FIRST_NEIGH";
				case StructureType::CUBIC_DIAMOND_SECOND_NEIGH: return "CUBIC_DIAMOND_SECOND_NEIGH";
				case StructureType::HEX_DIAMOND_FIRST_NEIGH: return "HEX_DIAMOND_FIRST_NEIGH";
				case StructureType::HEX_DIAMOND_SECOND_NEIGH: return "HEX_DIAMOND_SECOND_NEIGH";
                case StructureType::HEX_DIAMOND: return "HEX_DIAMOND";
                case StructureType::GRAPHENE: return "GRAPHENE";
                default: return "UNKNOWN";
            }
        } else {
            switch (structureType) {
                case static_cast<int>(CoordinationStructureType::COORD_OTHER): return "OTHER";
                case static_cast<int>(CoordinationStructureType::COORD_FCC): return "FCC";
                case static_cast<int>(CoordinationStructureType::COORD_HCP): return "HCP";
                case static_cast<int>(CoordinationStructureType::COORD_BCC): return "BCC";
                case static_cast<int>(CoordinationStructureType::COORD_CUBIC_DIAMOND): return "CUBIC_DIAMOND";
                case static_cast<int>(CoordinationStructureType::COORD_HEX_DIAMOND): return "HEX_DIAMOND";
                case static_cast<int>(StructureType::CUBIC_DIAMOND_FIRST_NEIGH): return "CUBIC_DIAMOND_FIRST_NEIGH";
				case static_cast<int>(StructureType::CUBIC_DIAMOND_SECOND_NEIGH): return "CUBIC_DIAMOND_SECOND_NEIGH";
				case static_cast<int>(StructureType::HEX_DIAMOND_FIRST_NEIGH): return "HEX_DIAMOND_FIRST_NEIGH";
				case static_cast<int>(StructureType::HEX_DIAMOND_SECOND_NEIGH): return "HEX_DIAMOND_SECOND_NEIGH";
                default: return "UNKNOWN";
            }
        }
    }

	std::string getAnalysisMethodName() const {
		switch(_identificationMode){
			case Mode::PTM: return "PTM";
			case Mode::CNA: return "CNA";
			case Mode::DIAMOND: return "DIAMOND";
			default: return "UNKNOWN";
		}
	}
    
    json buildMainListing() const{
		if(!_statisticsValid) calculateStructureStatistics();

		const int N = _context.atomCount();
		const double invN = (N > 0) ? (100.0 / static_cast<double>(N)) : 0.0;

		constexpr int K = static_cast<int>(StructureType::NUM_STRUCTURE_TYPES);
		std::vector<std::string> nameCache(K);
		std::vector<char> hasName(K, 0);

		auto getNameCached = [&](int st) -> const std::string& {
			const int idx = (0 <= st && st < K) ? st : static_cast<int>(StructureType::OTHER);
			if(!hasName[idx]){
				nameCache[idx] = getStructureTypeName(idx);
				hasName[idx] = 1;
			}
			return nameCache[idx];
		};

		int totalIdentified = 0;
		int unidentified = 0;
		auto itOther = _structureStatistics.find(static_cast<int>(StructureType::OTHER));
		if(itOther != _structureStatistics.end()){
			unidentified = itOther->second;
		}

		json mainListing = json::object();
		mainListing["total_atoms"] = N;
		mainListing["analysis_method"] = getAnalysisMethodName();

		for(const auto& [structureType, count] : _structureStatistics){
			std::string name = getNameCached(structureType);
			std::transform(name.begin(), name.end(), name.begin(), ::tolower);
			mainListing[name + "_count"] = count;
			mainListing[name + "_percentage"] = static_cast<double>(count) * invN;
			if(structureType != static_cast<int>(StructureType::OTHER) &&
				structureType != static_cast<int>(CoordinationStructureType::COORD_OTHER)){
				totalIdentified += count;
			}
		}

		mainListing["total_identified"] = totalIdentified;
		mainListing["total_unidentified"] = unidentified;
		mainListing["identification_rate"] = static_cast<double>(totalIdentified) * invN;
		mainListing["unique_structure_types"] = static_cast<int>(_structureStatistics.size());

		return mainListing;
	}



private:
	mutable std::map<int, int> _structureStatistics;
    mutable bool _statisticsValid = false;

	Mode _identificationMode;
	AnalysisContext& _context;
	bool _identifyPlanarDefects;
	std::mutex cluster_graph_mutex;
	
	float _rmsd;

	std::shared_ptr<ClusterGraph> _clusterGraph; 
	std::atomic<double> _maximumNeighborDistance;
};

}
