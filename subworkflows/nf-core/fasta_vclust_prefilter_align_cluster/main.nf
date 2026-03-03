include { VCLUST_PREFILTER } from '../../../modules/nf-core/vclust/prefilter/main'
include { VCLUST_ALIGN     } from '../../../modules/nf-core/vclust/align/main'
include { VCLUST_CLUSTER   } from '../../../modules/nf-core/vclust/cluster/main'

workflow FASTA_VCLUST_PREFILTER_ALIGN_CLUSTER {

    take:
    ch_fasta       // channel: [ val(meta), [ fasta ] ]
    save_alignment // boolean
    metric         // string
    tani           // float
    gani           // float
    ani            // float

    main:
    ch_versions = channel.empty()

    VCLUST_PREFILTER ( ch_fasta )
    ch_versions = ch_versions.mix(VCLUST_PREFILTER.out.versions.first())

    // Join to ensure in sync in case of multiple sequence files
    ch_align_in = VCLUST_PREFILTER.out.txt
        .join(ch_fasta)
        .multiMap { meta, txt, fasta ->
            fasta: [ meta, fasta ]
            filter: [ meta, txt ]
        }

    VCLUST_ALIGN ( ch_align_in.fasta, ch_align_in.filter, save_alignment )
    ch_versions = ch_versions.mix(VCLUST_ALIGN.out.versions.first())

    VCLUST_CLUSTER ( VCLUST_ALIGN.out.tsv, VCLUST_ALIGN.out.ids, metric, tani, gani, ani )
    ch_versions = ch_versions.mix(VCLUST_CLUSTER.out.versions.first())

    // Join to ensure in sync in case of multiple sequence files
    ch_out = VCLUST_CLUSTER.out.clusters
        .join(ch_fasta)
        .multiMap { meta, tsv, fasta ->
            fasta: [ meta, fasta ]
            tsv:   [ meta, tsv ]
        }

    emit:
    fasta     = ch_out.fasta           // channel: [ val(meta), [ fasta ] ]
    clusters  = ch_out.tsv             // channel: [ val(meta), [ tsv ] ]
    alignment = VCLUST_ALIGN.out.aln   // channel: [ val(meta), [ aln ] ]
    ids       = VCLUST_ALIGN.out.ids   // channel: [ val(meta), [ ids ] ]
    versions  = ch_versions            // channel: [ versions.yml ]
}
