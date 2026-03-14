//
// Run QC steps on BAM/CRAM files using Picard
//

include { PICARD_COLLECTMULTIPLEMETRICS } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTWGSMETRICS      } from '../../../modules/nf-core/picard/collectwgsmetrics/main'
include { PICARD_COLLECTHSMETRICS       } from '../../../modules/nf-core/picard/collecthsmetrics/main'

workflow BAM_QC_PICARD {
    take:
    ch_bam_bai_bait_target  // channel: [ val(meta), [bam], [bai], [bait_interval], [target_interval]]
    ch_fasta                // channel: [ val(meta), fasta ]
    ch_fasta_fai            // channel: [ val(meta), fasta_fai ]
    ch_fasta_dict           // channel: [ val(meta), fasta_dict ]
    ch_fasta_gzi            // channel: [ val(meta), fasta_gzi ]

    main:
    ch_coverage_metrics = channel.empty()

    ch_bam_bai = ch_bam_bai_bait_target.map{meta, bam, bai, _bait, _target -> return [meta,bam,bai]}

    PICARD_COLLECTMULTIPLEMETRICS( ch_bam_bai, ch_fasta, ch_fasta_fai )

    ch_bam_bai_bait_target_branched = ch_bam_bai_bait_target.branch {
        meta, bam, bai, bait, target ->
        hsmetrics  : bait != [] && target != []
            return [ meta, bam, bai, bait, target ]
        wgsmetrics : true
            return [ meta, bam, bai ]
    }

    PICARD_COLLECTHSMETRICS( ch_bam_bai_bait_target_branched.hsmetrics, ch_fasta, ch_fasta_fai, ch_fasta_dict, ch_fasta_gzi )
    ch_coverage_metrics = ch_coverage_metrics.mix(PICARD_COLLECTHSMETRICS.out.metrics)

    PICARD_COLLECTWGSMETRICS( ch_bam_bai_bait_target_branched.wgsmetrics, ch_fasta, ch_fasta_fai, [] )
    ch_coverage_metrics = ch_coverage_metrics.mix(PICARD_COLLECTWGSMETRICS.out.metrics)

    emit:
    coverage_metrics    = ch_coverage_metrics                       // channel: [ val(meta), [ coverage_metrics ] ]
    multiple_metrics    = PICARD_COLLECTMULTIPLEMETRICS.out.metrics // channel: [ val(meta), [ multiple_metrics ] ]
}
