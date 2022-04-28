//
// Run QC steps on BAM/CRAM files using Picard
//

params.options = [:]

include { PICARD_COLLECTMULTIPLEMETRICS } from '../../../modules/picardcollectmultiplemetrics/main'
include { PICARD_COLLECTWGSMETRICS      } from '../../../modules/picardcollectwgsmetrics/main'
include { PICARD_COLLECTHSMETRICS       } from '../../../modules/picardcollecthsmetrics/main'

workflow BAM_QC_PICARD {
    take:
    ch_bam_bai          // channel: [ val(meta), [ bam ], [bai/csi] ]
    ch_fasta            // channel: [ fasta ]
    ch_bait_interval    // channel: [ bait_interval ]
    ch_target_interval  // channel: [ target_interval ]

    main:
    ch_versions = Channel.empty()

    PICARD_COLLECTMULTIPLEMETRICS( ch_bam_bai, ch_fasta] )
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
    if (!ch_bait_interval.isEmpty() || !ch_target_interval.isEmpty()) {
        if (ch_bait_interval.isEmpty()) {
            throw new Error("Bait interval channel is empty")
        }
        if (ch_target_interval.isEmpty()) {
            throw new Error("Target interval channel is empty")
        }
        PICARD_COLLECTHSMETRICS( ch_bam_bai, ch_fasta, ch_bait_interval, ch_target_interval )
        ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions.first())
    } else {
        PICARD_COLLECTWGSMETRICS( ch_bam_bai, ch_fasta )
        ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions.first())
    }

    emit:
    hs_metrics          = PICARD_COLLECTHSMETRICS.out.hs_metrics    // channel: [ val(meta), [ hs_metrics ] ]
    wgs_metrics         = PICARD_COLLECTWGSMETRICS.out.metrics      // channel: [ val(meta), [ wgs_metrics ] ]
    multiple_metrics    = PICARD_COLLECTMULTIPLEMETRICS.out.metrics // channel: [ val(meta), [ multiple_metrics ] ]

    versions = ch_versions                                          // channel: [ versions.yml ]
}
