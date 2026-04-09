#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_analysis_context.h>

namespace Volt {

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
};

inline void normalizeReconstructedClusterGraphForExport(
    StructureAnalysis& analysis,
    AnalysisContext& context
){
    ReconstructedStateCanonicalizer::canonicalizeSymmetryPermutations(analysis, context);
}

}
