//
// Run QC steps on BAM/CRAM files using Picard
//

include { PICARD_COLLECTMULTIPLEMETRICS } from '../../../modules/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTWGSMETRICS      } from '../../../modules/picard/collectwgsmetrics/main'
include { PICARD_COLLECTHSMETRICS       } from '../../../modules/picard/collecthsmetrics/main'

workflow BAM_QC_PICARD {
    take:
    ch_bam              // channel: [ val(meta), [ bam ]]
    ch_fasta            // channel: [ fasta ]
    ch_fasta_fai        // channel: [ fasta_fai ]
    ch_bait_interval    // channel: [ bait_interval ]
    ch_target_interval  // channel: [ target_interval ]

    main:
    ch_versions = Channel.empty()
    ch_coverage_metrics = Channel.empty()

    PICARD_COLLECTMULTIPLEMETRICS( ch_bam, ch_fasta )
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())

    if (ch_bait_interval || ch_target_interval) {
        if (!ch_bait_interval) log.error("Bait interval channel is empty")
        if (!ch_target_interval) log.error("Target interval channel is empty")
        PICARD_COLLECTHSMETRICS( ch_bam, ch_fasta, ch_fasta_fai, ch_bait_interval, ch_target_interval )
        ch_coverage_metrics = ch_coverage_metrics.mix(PICARD_COLLECTHSMETRICS.out.metrics)
        ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions.first())
    } else {
        PICARD_COLLECTWGSMETRICS( ch_bam, ch_fasta )
        ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions.first())
        ch_coverage_metrics = ch_coverage_metrics.mix(PICARD_COLLECTWGSMETRICS.out.metrics)
    }

    emit:
    coverage_metrics    = ch_coverage_metrics                       // channel: [ val(meta), [ coverage_metrics ] ]
    multiple_metrics    = PICARD_COLLECTMULTIPLEMETRICS.out.metrics // channel: [ val(meta), [ multiple_metrics ] ]

    versions            = ch_versions                               // channel: [ versions.yml ]
}
