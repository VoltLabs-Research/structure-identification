#include <volt/analysis/structure_analysis.h>

namespace Volt {

void StructureAnalysis::identifyStructures(){
    if(usingPTM()){
        determineLocalStructuresWithPTM();
        computeMaximumNeighborDistanceFromPTM();
        return;
    }

    identifyStructuresCNA();
}

}
