#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_analysis_context.h>

#include <cstddef>
#include <vector>

namespace Volt {

struct CrystalTopologyEntry;

// How the neighbor shells of the crystalline atoms were put in export order: matched to their cluster's geometry,
// permuted by their symmetry, or kept in their local order when neither worked.
struct NeighborShellOrder{
    std::size_t geometryRemappedAtoms = 0;
    std::size_t permutationRemappedAtoms = 0;
    std::size_t localOrderAtoms = 0;
};

class ReconstructedStateCanonicalizer{
public:
    static void canonicalizeConnectedStructureSymmetries(
        StructureAnalysis& analysis,
        AnalysisContext& context
    );

    static void canonicalizeNeighborShellsToExportConvention(
        StructureAnalysis& analysis,
        AnalysisContext& context
    );

    static void canonicalizeSymmetryPermutations(
        StructureAnalysis& analysis,
        AnalysisContext& context
    );

private:
    static const CrystalTopologyEntry* topologyForStructure(
        const StructureAnalysis& analysis,
        int structureType
    );
    static int selectInitialSymmetry(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        int atomIndex,
        int structureType
    );
    static bool tryResolveNeighborSymmetry(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        int currentAtomIndex,
        int currentSymmetry,
        int neighborAtomIndex,
        int neighborIndex,
        int structureType,
        int& outNeighborSymmetry
    );
    static bool shouldVisitAtom(
        const AnalysisContext& context,
        std::size_t atomIndex,
        bool requireAssignedCluster
    );
    static bool shouldTraverseNeighbor(
        const AnalysisContext& context,
        int currentClusterId,
        int neighborAtomIndex,
        bool requireAssignedCluster
    );
    static double normalizedDotProduct(const Vector3& left, const Vector3& right);
    static bool solveMaximumWeightAssignment(
        const std::vector<std::vector<double>>& weights,
        std::vector<int>& assignment
    );
    static bool remapNeighborShellByClusterGeometry(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        std::size_t atomIndex,
        int structureType,
        int start,
        int exportableCount,
        const int* oldIndices,
        const CrystalTopologyEntry& topology,
        std::vector<int>& localSlotToExportSlot
    );
    static void canonicalizeSymmetryPermutationsImpl(
        StructureAnalysis& analysis,
        AnalysisContext& context,
        bool requireAssignedCluster
    );
};

inline void normalizeReconstructedClusterGraphForExport(
    StructureAnalysis& analysis,
    AnalysisContext& context
){
    ReconstructedStateCanonicalizer::canonicalizeSymmetryPermutations(analysis, context);
}

}
