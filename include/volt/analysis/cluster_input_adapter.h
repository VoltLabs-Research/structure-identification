#pragma once

#include <volt/analysis/analysis_context.h>
#include <volt/analysis/structure_analysis.h>

namespace Volt{

class ClusterInputAdapter{
public:
    virtual ~ClusterInputAdapter() = default;

    virtual void prepare(StructureAnalysis& analysis, AnalysisContext& context) = 0;
};

}
