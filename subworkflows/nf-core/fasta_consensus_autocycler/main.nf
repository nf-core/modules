
include { AUTOCYCLER_COMPRESS  } from '../../../modules/nf-core/autocycler/compress/main'
include { AUTOCYCLER_CLUSTER   } from '../../../modules/nf-core/autocycler/cluster/main'
include { AUTOCYCLER_TRIM      } from '../../../modules/nf-core/autocycler/trim/main'
include { AUTOCYCLER_RESOLVE   } from '../../../modules/nf-core/autocycler/resolve/main'
include { AUTOCYCLER_COMBINE   } from '../../../modules/nf-core/autocycler/combine/main'

workflow FASTA_CONSENSUS_AUTOCYCLER {

    take:
    ch_grouped_contigs // channel: [ val(meta), [ fasta, fasta, ... ] ]

    main:

    AUTOCYCLER_COMPRESS (
        ch_grouped_contigs
    )

    AUTOCYCLER_CLUSTER (
        AUTOCYCLER_COMPRESS.out.gfa
    )

    AUTOCYCLER_CLUSTER.out.clusters
        .flatMap { meta, gfa_list ->
            ( gfa_list instanceof List ? gfa_list: [gfa_list] )
            .collect { gfa -> [meta, gfa] } // Separate gfas, each with a meta
        }
        .map { meta, file ->
            def cluster_id = (file.parent.name) // Add cluster id to meta
            [ meta + [cluster: cluster_id], file ]
        }
        .set{ ch_clusters } // channel: [ val(meta), gfa]

    AUTOCYCLER_TRIM(
        ch_clusters
    )

    AUTOCYCLER_RESOLVE(
        AUTOCYCLER_TRIM.out.gfa
    )

    // Group clusters based on meta, ignoring cluster id
    AUTOCYCLER_RESOLVE.out.resolved
        .map{   meta, file ->
            def new_meta = meta.clone()
            new_meta.remove("cluster")
            tuple( new_meta, file)
        }
        .groupTuple()
        .set{ ch_assembly_clusters } // channel: [ meta, [ cluster1, cluster2, ... ] ]

    AUTOCYCLER_COMBINE(
        ch_assembly_clusters
    )

    ch_consensus_assembly       = AUTOCYCLER_COMBINE.out.fasta // channel: [ val(meta), fasta ]
    ch_consensus_assembly_graph = AUTOCYCLER_COMBINE.out.gfa   // channel: [ val(meta), gfa ]


    emit:
    consensus_assembly       = ch_consensus_assembly       // channel: [ val(meta), fasta ]
    consensus_assembly_graph = ch_consensus_assembly_graph // channel: [ val(meta), gfa ]
}
