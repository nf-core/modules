//
// Optionally concatenate FASTQ files, align with minimap2, index BAM, and run BAM stats.
//

include { CAT_FASTQ          } from '../../../modules/nf-core/cat/fastq/main'
include { MINIMAP2_ALIGN     } from '../../../modules/nf-core/minimap2/align/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'

workflow FASTQ_ALIGN_MINIMAP2 {

    take:
    ch_reads
    // [ val(meta), path(reads) ]
    // reads may be one FASTQ or a list of FASTQs. Use meta.single_end = true for ONT.

    ch_fasta_fai
    // [ val(meta2), path(fasta), path(fai) ]
    // reference can be FASTA or prebuilt minimap2 .mmi

    val_bam_index_extension
    // "bai" or "csi"

    val_cigar_bam
    // true/false; passed to minimap2/align as long-CIGAR BAM handling

    val_skip_cat
    // true: pass input FASTQ(s) directly to minimap2
    // false: concatenate with CAT_FASTQ before minimap2

    main:

    ch_fasta = ch_fasta_fai.map { meta, fasta, _fai -> [ meta, fasta ] }

    ch_merged_fastq = channel.empty()

    if (val_skip_cat) {
        ch_reads_for_alignment = ch_reads
    } else {
        CAT_FASTQ(ch_reads)
        ch_reads_for_alignment = CAT_FASTQ.out.reads
        ch_merged_fastq        = CAT_FASTQ.out.reads
    }

    MINIMAP2_ALIGN(
        ch_reads_for_alignment,
        ch_fasta,
        true,
        val_bam_index_extension,
        false,
        val_cigar_bam
    )

    MINIMAP2_ALIGN.out.bam
        .join(MINIMAP2_ALIGN.out.index, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .set { ch_bam_bai }

    BAM_STATS_SAMTOOLS(ch_bam_bai, ch_fasta_fai)

    emit:
    merged_fastq = ch_merged_fastq
    bam          = MINIMAP2_ALIGN.out.bam
    index        = MINIMAP2_ALIGN.out.index
    stats        = BAM_STATS_SAMTOOLS.out.stats
    flagstat     = BAM_STATS_SAMTOOLS.out.flagstat
    idxstats     = BAM_STATS_SAMTOOLS.out.idxstats
}
