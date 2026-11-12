// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

/*
 * TAPS methylation conversion subworkflow
 *
 * Uses rastair to assess C->T conversion as a readout for methylation in a genome-wide basis.
 * By default runs rastair 0.8 (mbias -> mbiasparser -> call -> methylkit). With `use_rastair2`,
 * runs rastair version 2 instead (mbias -> call), which writes a VCF (methylation, SNVs) and a CpG BED, plus
 * bcftools stats. rastair2 mbias needs paired-end reads, so single-end samples are called without trimming.
 * With `trims`, both callers use those bases to ignore and the M-bias steps are skipped.
 * Note: rastair version 2 is free for non-commercial use only (University of Oxford licence).
 */

include { RASTAIR_MBIAS             } from '../../../modules/nf-core/rastair/mbias/main'
include { RASTAIR_MBIASPARSER       } from '../../../modules/nf-core/rastair/mbiasparser/main'
include { RASTAIR_CALL              } from '../../../modules/nf-core/rastair/call/main'
include { RASTAIR_METHYLKIT         } from '../../../modules/nf-core/rastair/methylkit/main'
include { RASTAIR2_MBIAS            } from '../../../modules/nf-core/rastair2/mbias/main'
include { RASTAIR2_CALL             } from '../../../modules/nf-core/rastair2/call/main'
include { BCFTOOLS_STATS            } from '../../../modules/nf-core/bcftools/stats/main'

workflow BAM_TAPS_CONVERSION {

    take:
    ch_bam                 // channel: [ val(meta), [ bam ] ] ## BAM from alignment
    ch_bai                 // channel: [ val(meta), [ bai ] ] ## BAI from alignment
    ch_fasta               // channel: [ val(meta), [ fa ] ]
    ch_fasta_index         // channel: [ val(meta), [ fa.fai ] ]
    use_rastair2           // boolean: run rastair version 2 instead of rastair 0.8
    trims                  // list: [ nOT, nOB ] bases to ignore, each 'r1_start,r1_end,r2_start,r2_end'; [] to use the M-bias estimate

    main:
    ch_rastair_mbias = channel.empty()
    ch_rastair_call  = channel.empty()
    ch_methylkit     = channel.empty()
    ch_vcf           = channel.empty()
    ch_vcf_tbi       = channel.empty()
    ch_bed           = channel.empty()
    ch_bed_tbi       = channel.empty()
    ch_bcftools_stats = channel.empty()
    ch_mbias_html    = channel.empty()
    ch_mbias_cutoffs = channel.empty()

    if (trims && !(trims.size() == 2 && trims.every { trim -> trim ==~ /\d+,\d+,\d+,\d+/ })) {
        error("bam_taps_conversion: trims must be [] or [ nOT, nOB ], each as 'r1_start,r1_end,r2_start,r2_end' (got ${trims})")
    }

    if (use_rastair2) {
        // One reference for every sample, whether it is passed once or once per sample
        ch_reference = ch_fasta
            .combine(ch_fasta_index)
            .map { meta, fasta, _meta2, fai -> [ meta, fasta, fai ] }
            .first()

        ch_bam_bai = ch_bam.join(ch_bai)

        if (trims) {
            ch_trims_file = channel.of("${trims[0]}\n${trims[1]}\n")
                .collectFile(name: 'trims.cutoffs.txt')
                .first()
            ch_call_cutoffs = ch_bam_bai
                .combine(ch_trims_file)
                .map { meta, _bam, _bai, cutoffs -> [ meta, cutoffs ] }
        }
        else {
            RASTAIR2_MBIAS (
                ch_bam_bai.filter { meta, _bam, _bai -> !meta.single_end },
                ch_reference,
            )
            ch_mbias_html = RASTAIR2_MBIAS.out.html // channel: [ val(meta), html ]

            // Without --region, mbias writes cutoffs for every contig with enough data: use the longest contig
            ch_mbias_cutoffs = RASTAIR2_MBIAS.out.cutoffs
                .combine(ch_reference)
                .map { meta, cutoffs, _meta2, _fasta, fai ->
                    def files = cutoffs instanceof List ? cutoffs : [ cutoffs ]
                    def contig_lengths = fai.readLines().collectEntries { line ->
                        def fields = line.split('\t')
                        [ (fields[0]): fields[1] as long ]
                    }
                    def longest = files.max { cutoffs_file ->
                        def contig = contig_lengths.keySet()
                            .findAll { name -> cutoffs_file.name.endsWith(".${name}.cutoffs.txt") }
                            .max { name -> name.length() }
                        contig ? contig_lengths[contig] : 0L
                    }
                    [ meta, longest ]
                } // channel: [ val(meta), cutoffs.txt ]
            ch_call_cutoffs = ch_mbias_cutoffs
        }

        RASTAIR2_CALL (
            ch_bam_bai
                .join(ch_call_cutoffs, remainder: true)
                .map { meta, bam, bai, cutoffs -> [ meta, bam, bai, cutoffs ?: [] ] },
            ch_reference,
        )
        ch_vcf     = RASTAIR2_CALL.out.vcf     // channel: [ val(meta), vcf.gz ]
        ch_vcf_tbi = RASTAIR2_CALL.out.tbi     // channel: [ val(meta), vcf.gz.tbi ]
        ch_bed     = RASTAIR2_CALL.out.bed     // channel: [ val(meta), bed.gz ]
        ch_bed_tbi = RASTAIR2_CALL.out.bed_tbi // channel: [ val(meta), bed.gz.tbi ]

        BCFTOOLS_STATS (
            ch_vcf.join(ch_vcf_tbi),
            [ [:], [] ],
            [ [:], [] ],
            [ [:], [] ],
            [ [:], [] ],
            [ [:], [] ],
        )
        ch_bcftools_stats = BCFTOOLS_STATS.out.stats // channel: [ val(meta), stats.txt ]
    }
    else {
        if (trims) {
            ch_trims = ch_bam.map { meta, _bam -> [ meta, trims[0], trims[1] ] }
        }
        else {
            RASTAIR_MBIAS (
                ch_bam,
                ch_bai,
                ch_fasta,
                ch_fasta_index,
            )
            ch_rastair_mbias = RASTAIR_MBIAS.out.txt // channel: [ val(meta), txt ]

            RASTAIR_MBIASPARSER (
                ch_rastair_mbias
            )
            ch_trims = RASTAIR_MBIASPARSER.out.mbias_processed_str // channel: [ val(meta), nOT_clip, nOB_clip ]
        }

        RASTAIR_CALL (
            ch_bam,
            ch_bai,
            ch_fasta,
            ch_fasta_index,
            ch_trims.map{ meta, nOT_clip, _nOB_clip -> [ meta, nOT_clip ] },
            ch_trims.map{ meta, _nOT_clip, nOB_clip -> [ meta, nOB_clip ] },
        )
        ch_rastair_call = RASTAIR_CALL.out.txt // channel: [ val(meta), txt ]

        RASTAIR_METHYLKIT (
            ch_rastair_call
        )
        ch_methylkit = RASTAIR_METHYLKIT.out.methylkit // channel
    }

    emit:
    mbias          = ch_rastair_mbias         // channel: [ val(meta), path("*.txt") ]    rastair 0.8
    call           = ch_rastair_call          // channel: [ val(meta), path("*.txt") ]    rastair 0.8
    methylkit      = ch_methylkit             // channel: [ val(meta), path("*.txt.gz") ] rastair 0.8
    vcf            = ch_vcf                   // channel: [ val(meta), path("*.vcf.gz") ] rastair 2
    vcf_tbi        = ch_vcf_tbi               // channel: [ val(meta), path("*.vcf.gz.tbi") ] rastair 2
    bed            = ch_bed                   // channel: [ val(meta), path("*.bed.gz") ] rastair 2
    bed_tbi        = ch_bed_tbi               // channel: [ val(meta), path("*.bed.gz.tbi") ] rastair 2
    bcftools_stats = ch_bcftools_stats        // channel: [ val(meta), path("*stats.txt") ] rastair 2
    mbias_html     = ch_mbias_html            // channel: [ val(meta), path("*.mbias.html") ] rastair 2
    mbias_cutoffs  = ch_mbias_cutoffs         // channel: [ val(meta), path("*.cutoffs.txt") ] rastair 2
}
