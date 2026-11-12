// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

//
// Align reads with bwa-mem3 (plain, bisulfite/EM-seq `--meth` or TAPS `--meth=taps`, set via BWAMEM3_MEM ext.args),
// coordinate-sort and index, optionally mark duplicates, and collect samtools stats
//

include { BWAMEM3_MEM           } from '../../../modules/nf-core/bwamem3/mem/main'
include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from '../../nf-core/bam_stats_samtools/main'

workflow FASTQ_ALIGN_DEDUP_BWAMEM3 {
    take:
    ch_reads // channel: [ val(meta), [ reads ] ]
    ch_fasta_fai // channel: [ val(meta), fasta, fai ]
    ch_bwamem3_index // channel: [ val(meta), [ bwa-mem3 index ] ]
    skip_deduplication // boolean: whether to skip marking duplicates

    main:
    ch_picard_metrics = channel.empty()

    // Align, coordinate-sort and write a BAI in one step
    BWAMEM3_MEM(
        ch_reads,
        ch_bwamem3_index,
        ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] },
        true,
        'bai',
    )
    ch_alignment = BWAMEM3_MEM.out.aligned
    ch_alignment_index = BWAMEM3_MEM.out.index

    if (!skip_deduplication) {
        PICARD_MARKDUPLICATES(
            ch_alignment,
            ch_fasta_fai,
        )
        SAMTOOLS_INDEX(PICARD_MARKDUPLICATES.out.bam)
        ch_alignment = PICARD_MARKDUPLICATES.out.bam
        ch_alignment_index = SAMTOOLS_INDEX.out.index
        ch_picard_metrics = PICARD_MARKDUPLICATES.out.metrics
    }

    // samtools stats, flagstat and idxstats on the final BAM
    BAM_STATS_SAMTOOLS(
        ch_alignment.join(ch_alignment_index),
        ch_fasta_fai,
    )

    ch_multiqc_files = ch_picard_metrics
        .mix(BAM_STATS_SAMTOOLS.out.stats)
        .mix(BAM_STATS_SAMTOOLS.out.flagstat)
        .mix(BAM_STATS_SAMTOOLS.out.idxstats)
        .collect { _meta, file -> file }

    emit:
    bam               = ch_alignment // channel: [ val(meta), bam ]
    index             = ch_alignment_index // channel: [ val(meta), bai ]
    samtools_flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), flagstat ]
    samtools_stats    = BAM_STATS_SAMTOOLS.out.stats // channel: [ val(meta), stats ]
    samtools_idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), idxstats ]
    picard_metrics    = ch_picard_metrics // channel: [ val(meta), metrics ]
    multiqc           = ch_multiqc_files // channel: [ files ]
}
