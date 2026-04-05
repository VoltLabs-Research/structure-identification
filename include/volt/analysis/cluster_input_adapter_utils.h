#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/symmetry_utils.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>

namespace Volt::ClusterInputAdapterUtils{

template<typename ShouldAssignMask>
void prepareSymmetryAwareClusterInputs(
    StructureAnalysis& analysis,
    AnalysisContext& context,
    bool resetAllowedSymmetryMasks,
    ShouldAssignMask&& shouldAssignMask
){
    analysis.setClusterRuleProvider(nullptr);

    if(!context.atomAllowedSymmetryMasks){
        context.atomAllowedSymmetryMasks = std::make_shared<ParticleProperty>(
            context.atomCount(),
            DataType::Int64,
            1,
            0,
            true
        );
    }

    std::fill(
        context.atomSymmetryPermutations->dataInt(),
        context.atomSymmetryPermutations->dataInt() + context.atomSymmetryPermutations->size(),
        -1
    );

    if(resetAllowedSymmetryMasks){
        std::fill(
            context.atomAllowedSymmetryMasks->dataInt64(),
            context.atomAllowedSymmetryMasks->dataInt64() + context.atomAllowedSymmetryMasks->size(),
            0
        );
    }

    for(std::size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
        const int structureType = context.structureTypes->getInt(atomIndex);
        if(!shouldAssignMask(atomIndex, structureType)){
            continue;
        }

        context.atomAllowedSymmetryMasks->setInt64(
            atomIndex,
            static_cast<std::int64_t>(
                AnalysisSymmetryUtils::fullSymmetryMask(analysis.symmetryPermutationCount(structureType))
            )
        );
    }
}

}
