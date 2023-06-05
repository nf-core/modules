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

    main:
    ch_versions = Channel.empty()
    ch_coverage_metrics = Channel.empty()

    ch_bam_bai = ch_bam_bai_bait_target.map{meta, bam, bai, bait, target -> return [meta,bam,bai]}

    PICARD_COLLECTMULTIPLEMETRICS( ch_bam_bai, ch_fasta, ch_fasta_fai )
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())

    ch_bam_bai_bait_target_branched = ch_bam_bai_bait_target.branch {
        hsmetrics  : it.size == 5 && it[3] != [] && it[4] != []
            return it
        wgsmetrics : true
            return [ it[0], it[1], it[2] ]
    }

    PICARD_COLLECTHSMETRICS( ch_bam_bai_bait_target_branched.hsmetrics, ch_fasta, ch_fasta_fai, ch_fasta_dict )
    ch_coverage_metrics = ch_coverage_metrics.mix(PICARD_COLLECTHSMETRICS.out.metrics)
    ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions.first())

    PICARD_COLLECTWGSMETRICS( ch_bam_bai_bait_target_branched.wgsmetrics, ch_fasta, ch_fasta_fai, [] )
    ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions.first())
    ch_coverage_metrics = ch_coverage_metrics.mix(PICARD_COLLECTWGSMETRICS.out.metrics)

    emit:
    coverage_metrics    = ch_coverage_metrics                       // channel: [ val(meta), [ coverage_metrics ] ]
    multiple_metrics    = PICARD_COLLECTMULTIPLEMETRICS.out.metrics // channel: [ val(meta), [ multiple_metrics ] ]

    versions            = ch_versions                               // channel: [ versions.yml ]
}
