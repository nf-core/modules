//
// Run QC steps on BAM/CRAM files using Picard
//

include { PICARD_COLLECTMULTIPLEMETRICS } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTWGSMETRICS      } from '../../../modules/nf-core/picard/collectwgsmetrics/main'
include { PICARD_COLLECTHSMETRICS       } from '../../../modules/nf-core/picard/collecthsmetrics/main'

workflow BAM_QC_PICARD {
    take:
    ch_bam_bai          // channel: [mandatory] [ val(meta), [ path(bam) ], [ path(bai) ] ]
    ch_fasta            // channel: [mandatory] [ val(meta2), path(fasta) ]
    ch_fasta_fai        // channel: [optional]  [ val(meta2), path(fasta_fai) ]
    ch_bait_interval    // channel: [optional]  [ path(bait_interval) ]
    ch_target_interval  // channel: [optional]  [ path(target_interval) ]
    ch_interval_list    // channel: [optional]  [ path(interval_list) ]

    main:
    ch_versions = Channel.empty()
    ch_coverage_metrics = Channel.empty()

    PICARD_COLLECTMULTIPLEMETRICS( ch_bam_bai, ch_fasta, ch_fasta_fai )
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())

    if (ch_bait_interval || ch_target_interval) {
        if (!ch_bait_interval) log.error("Bait interval channel is empty")
        if (!ch_target_interval) log.error("Target interval channel is empty")
        PICARD_COLLECTHSMETRICS( ch_bam_bai, ch_fasta, ch_fasta_fai, ch_bait_interval, ch_target_interval )
        ch_coverage_metrics = ch_coverage_metrics.mix(PICARD_COLLECTHSMETRICS.out.metrics)
        ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions.first())
    } else {
        PICARD_COLLECTWGSMETRICS( ch_bam_bai, ch_fasta, ch_fasta_fai, ch_interval_list )
        ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions.first())
        ch_coverage_metrics = ch_coverage_metrics.mix(PICARD_COLLECTWGSMETRICS.out.metrics)
    }

    emit:
    coverage_metrics     = ch_coverage_metrics                       // channel: [ val(meta), [ path(coverage_metrics) ] ]
    multiple_metrics     = PICARD_COLLECTMULTIPLEMETRICS.out.metrics // channel: [ val(meta), [ path(multiple_metrics) ] ]
    multiple_metrics_pdf = PICARD_COLLECTMULTIPLEMETRICS.out.pdf     // channel: [ val(meta), [ path(multiple_metrics_pdf) ] ]

    versions             = ch_versions                               // channel: [ path(versions.yml) ]
}
