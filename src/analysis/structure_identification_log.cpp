#include <volt/analysis/structure_identification_log.h>

#include <volt/analysis/cluster_graph_builder.h>
#include <volt/analysis/cluster_graph_io.h>
#include <volt/analysis/reconstructed_analysis_pipeline.h>
#include <volt/analysis/reconstructed_state_canonicalizer.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/core/reconstructed_structure.h>
#include <volt/structures/crystal_structure_types.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cstdint>
#include <vector>

namespace Volt{

void StructureIdentificationLog::dump(const ClusterGrowth& growth){
    spdlog::info("  Grew {} clusters covering {} atoms in {}", growth.clusters, growth.atoms, phaseElapsed());
    const std::size_t unclusteredSeeds = growth.crystallineAtoms - growth.atoms;
    if(unclusteredSeeds > 0){
        spdlog::info("  {} atoms have no allowed symmetry permutation; left unclustered as Other", unclusteredSeeds);
    }
}

void StructureIdentificationLog::dump(const ClusterConnections& connections){
    const TrimmedNeighbors& trimmed = connections.trimmedNeighbors;
    if(trimmed.neighbors > 0){
        spdlog::warn(
            "  Trimmed {} appended neighbors across {} atoms to respect MAX_NEIGHBORS={}",
            trimmed.neighbors, trimmed.atoms, static_cast<int>(MAX_NEIGHBORS)
        );
    }
    spdlog::info("  Connected clusters: {} transitions in {}", connections.transitions, phaseElapsed());
}

void StructureIdentificationLog::dump(const SuperClusterBridges& bridges){
    if(bridges.transitions > 0){
        spdlog::info("  Bridged {} cluster pairs across defect clusters", bridges.transitions);
    }
}

void StructureIdentificationLog::dump(const DissolvedClusters& dissolved){
    spdlog::info("  Marked {} clusters smaller than {} atoms as Other", dissolved.clusters, dissolved.minClusterSize);
}

void StructureIdentificationLog::dump(const NeighborShellOrder& order){
    spdlog::info(
        "  Neighbor shells in export order: {} atoms by cluster geometry, {} by symmetry permutation, {} kept in local order",
        order.geometryRemappedAtoms, order.permutationRemappedAtoms, order.localOrderAtoms
    );
}

void StructureIdentificationLog::dump(const AnnotatedDump& annotatedDump){
    spdlog::info("  Annotated dump written to {} in {}", annotatedDump.path, phaseElapsed());
}

void StructureIdentificationLog::dump(const ClusterGraphTables& tables){
    spdlog::info("  Cluster graph written: {} clusters, {} transitions", tables.clusters, tables.transitions);
}

void StructureIdentificationLog::dump(const NeighborLatticeTable& neighborLattice){
    spdlog::info(
        "  Neighbor lattice written to {} ({} atoms) in {}",
        neighborLattice.path, neighborLattice.atoms, phaseElapsed()
    );
}

void StructureIdentificationLog::dump(const ReconstructedStructureContext& context){
    const std::size_t atomCount = context.atomCount();
    const std::size_t neighborEntries = context.neighborIndices ? context.neighborIndices->size() : 0;
    spdlog::info(
        "  Neighbor lattice loaded: {} atoms, {:.1f} neighbors per atom, max neighbor distance {:.3f} in {}",
        atomCount,
        atomCount > 0 ? static_cast<double>(neighborEntries) / static_cast<double>(atomCount) : 0.0,
        context.maximumNeighborDistance,
        phaseElapsed()
    );
}

void StructureIdentificationLog::dump(const ClusterGraphImport& graphImport){
    spdlog::info(
        "  Cluster graph imported: {} clusters, {} transitions in {}",
        graphImport.graph.clusters().size() - graphImport.clustersBefore,
        graphImport.graph.clusterTransitions().size() - graphImport.transitionsBefore,
        phaseElapsed()
    );
}

void StructureIdentificationLog::dump(const BuiltClusters& built){
    const auto& clusters = built.graph.clusters();
    // The cluster graph always holds the null cluster 0.
    const std::size_t clusterCount = clusters.empty() ? 0 : clusters.size() - 1;
    if(built.latticeSites){
        const std::vector<int>& sites = *built.latticeSites;
        const auto clusteredSites = std::ranges::count_if(sites, [&](int atom){
            return built.atomClusters.getInt(atom) != 0;
        });
        spdlog::info("Built {} clusters in {} ({} of {} lattice sites in clusters, {:.1f}%)", clusterCount, phaseElapsed(),
            clusteredSites, sites.size(), percentOf(static_cast<std::uint64_t>(clusteredSites), sites.size()));
        return;
    }
    if(!built.smallClustersDissolved){
        spdlog::info("Built {} clusters in {}", clusterCount, phaseElapsed());
        return;
    }
    const auto dissolvedClusters = std::ranges::count_if(clusters, [](const Cluster* cluster){
        return cluster && cluster->id != 0 && cluster->structure == LATTICE_OTHER;
    });
    spdlog::info(
        "Built {} clusters in {} ({} with fewer than {} atoms dissolved)",
        clusterCount, phaseElapsed(), dissolvedClusters, ClusterBuilder::kDefaultSmallClusterSize
    );
}

void StructureIdentificationLog::dump(const AnalysisPipelineUtils::ClusterOutputs&){
    spdlog::info("Cluster outputs written in {}", phaseElapsed());
}

}
