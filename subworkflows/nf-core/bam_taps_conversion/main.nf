/*
 * TAPS methylation conversion subworkflow
 *
 * Uses Rastair to assess C->T conversion as a readout for methylation in a genome-wide basis
 */

include { RASTAIR_MBIAS             } from '../../../modules/nf-core/rastair/mbias/main'
include { RASTAIR_MBIASPARSER       } from '../../../modules/nf-core/rastair/mbiasparser/main'
include { RASTAIR_CALL              } from '../../../modules/nf-core/rastair/call/main'
include { RASTAIR_METHYLKIT         } from '../../../modules/nf-core/rastair/methylkit/main'

workflow BAM_TAPS_CONVERSION {

    take:
    ch_bam                 // channel: [ val(meta), [ bam ] ] ## BAM from alignment
    ch_bai                 // channel: [ val(meta), [ bai ] ] ## BAI from alignment
    ch_fasta               // channel: [ val(meta), [ fa ] ]
    ch_fasta_index         // channel: [ val(meta), [ fa.fai ] ]

    main:
    ch_rastair_mbias = channel.empty()
    ch_rastair_call  = channel.empty()
    ch_versions      = channel.empty()

    log.info "Running TAPS conversion module with Rastair to assess C->T conversion as a readout for methylation."

    RASTAIR_MBIAS (
        ch_bam,
        ch_bai,
        ch_fasta,
        ch_fasta_index,
    )
    ch_rastair_mbias = RASTAIR_MBIAS.out.txt // channel: [ val(meta), txt ]
    ch_versions      = ch_versions.mix(RASTAIR_MBIAS.out.versions)

    RASTAIR_MBIASPARSER (
        ch_rastair_mbias
    )
    ch_rastair_mbiasparser = RASTAIR_MBIASPARSER.out.mbias_processed_str // channel: [ val(meta), nOT_clip, nOB_clip ]
    ch_versions             = ch_versions.mix(RASTAIR_MBIASPARSER.out.versions)

    RASTAIR_CALL (
        ch_bam,
        ch_bai,
        ch_fasta,
        ch_fasta_index,
        ch_rastair_mbiasparser.map{ meta, nOT_clip, _nOB_clip -> [ meta, nOT_clip ] },
        ch_rastair_mbiasparser.map{ meta, _nOT_clip, nOB_clip -> [ meta, nOB_clip ] },
    )
    ch_rastair_call = RASTAIR_CALL.out.txt // channel: [ val(meta), txt ]
    ch_versions     = ch_versions.mix(RASTAIR_CALL.out.versions)

    RASTAIR_METHYLKIT (
        ch_rastair_call
    )
    ch_methylkit = RASTAIR_METHYLKIT.out.methylkit // channel
    ch_versions  = ch_versions.mix(RASTAIR_METHYLKIT.out.versions)

    emit:
    mbias        = ch_rastair_mbias         // channel: [ val(meta), path("*.txt") ]
    call         = ch_rastair_call          // channel: [ val(meta), path("*.txt") ]
    methylkit    = ch_methylkit             // channel: [ val(meta), path("*.txt.gz") ]
    versions     = ch_versions              // channel: path("*.version.txt")
}
