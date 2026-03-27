#pragma once

#include <volt/core/volt.h>
#include <volt/core/particle_property.h>
#include <volt/core/simulation_cell.h>
#include <volt/structures/cluster_graph.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/analysis/analysis_context.h>
#include <volt/core/lammps_parser.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace Volt{

class StructureAnalysisCrystalInfo{
public:
	virtual ~StructureAnalysisCrystalInfo() = default;

	virtual int findClosestSymmetryPermutation(int structureType, const Matrix3& rotation) const = 0;
	virtual int coordinationNumber(int structureType) const = 0;
	virtual int commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const = 0;
	virtual int symmetryPermutationCount(int structureType) const = 0;
	virtual int symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const = 0;
	virtual const Matrix3& symmetryTransformation(int structureType, int symmetryIndex) const = 0;
	virtual int symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const = 0;
	virtual const Vector3& latticeVector(int structureType, int latticeVectorIndex) const = 0;
};

struct PTMComputationData{
	std::shared_ptr<ParticleProperty> rmsd;
	std::shared_ptr<ParticleProperty> orientations;
	std::shared_ptr<ParticleProperty> correspondences;
};

class StructureAnalysis{
public:
	StructureAnalysis(StructureContext& context);

	void appendNeighbors(const std::vector<std::vector<int>>& extras);

	StructureContext& context(){
		return _context;
	}

	const StructureContext& context() const{
		return _context;
	}

	double maximumNeighborDistance() const{
		return _context.maximumNeighborDistance;
	}

	int numberOfNeighbors(int atomIndex) const {
		assert(_context.neighborCounts);
		return _context.neighborCounts->getInt(atomIndex);
	}
	
	int getNeighbor(int centralAtomIndex, int neighborListIndex) const;
	int findNeighbor(int centralAtomIndex, int neighborAtomIndex) const;

	const ClusterGraph& clusterGraph() const{
		return *_clusterGraph;
	}

	ClusterGraph& clusterGraph(){
		return *_clusterGraph;
	}

	// PTM
	void computeMaximumNeighborDistanceFromPTM();
	PTMComputationData determineLocalStructuresWithPTM(double rmsdCutoff);

	Cluster* atomCluster(int atomIndex) const{
		return clusterGraph().findCluster(_context.atomClusters->getInt(atomIndex));
	}
	
	int findClosestSymmetryPermutation(int structureType, const Matrix3& rotation);
	int coordinationNumber(int structureType) const;
	int commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const;
	int symmetryPermutationCount(int structureType) const;
	int symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const;
	const Matrix3& symmetryTransformation(int structureType, int symmetryIndex) const;
	int symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const;

	void identifyStructuresCNA();

	const Vector3& latticeVector(int structureType, int latticeVectorIndex) const;

	// Returns the ideal lattice vector associated with a neighbor bond
	const Vector3& neighborLatticeVector(int centralAtomIndex, int neighborIndex) const;

	void setNeighborLatticeVectorOverrides(std::vector<Vector3> overrides, std::size_t stride);
	void setCrystalInfoProvider(std::shared_ptr<const StructureAnalysisCrystalInfo> crystalInfoProvider);
	bool hasCrystalInfoProvider() const{
		return static_cast<bool>(_crystalInfoProvider);
	}
	bool hasNeighborLatticeVectorOverrides() const{
		return !_neighborLatticeVectorOverrides.empty();
	}
	const std::vector<Vector3>& neighborLatticeVectorOverrides() const{
		return _neighborLatticeVectorOverrides;
	}
	std::size_t neighborLatticeVectorOverrideStride() const{
		return _neighborLatticeVectorOverrideStride;
	}

private:
	AnalysisContext& analysisContext();
	const AnalysisContext& analysisContext() const;

	StructureContext& _context;
	std::shared_ptr<ClusterGraph> _clusterGraph; 
	std::shared_ptr<const StructureAnalysisCrystalInfo> _crystalInfoProvider;
	std::vector<Vector3> _neighborLatticeVectorOverrides;
	std::size_t _neighborLatticeVectorOverrideStride = 0;
};

}
