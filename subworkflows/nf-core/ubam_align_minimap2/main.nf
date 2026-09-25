//
// Optionally concatenate uBAM files, align with minimap2, index BAM, and run BAM stats.
//

include { SAMTOOLS_CAT      } from '../../../modules/nf-core/samtools/cat/main'
include { MINIMAP2_ALIGN    } from '../../../modules/nf-core/minimap2/align/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'

workflow UBAM_ALIGN_MINIMAP2 {

    take:
    ch_ubam
    // [ val(meta), path(ubam) ]
    // ubam may be one BAM/uBAM file or a list of BAM/uBAM files.

    ch_fasta_fai
    // [ val(meta2), path(fasta), path(fai) ]

    val_bam_index_extension
    // "bai" or "csi"

    val_cigar_bam
    // true/false

    val_skip_cat
    // true: pass one uBAM directly to minimap2
    // false: concatenate uBAM files with samtools/cat before minimap2

    main:

    ch_fasta = ch_fasta_fai.map { meta, fasta, _fai -> [ meta, fasta ] }

    ch_cat_ubam = channel.empty()

    if (val_skip_cat) {
        ch_ubam_for_alignment = ch_ubam
    } else {
        SAMTOOLS_CAT(ch_ubam)
        ch_ubam_for_alignment = SAMTOOLS_CAT.out.bam
        ch_cat_ubam           = SAMTOOLS_CAT.out.bam
    }

    MINIMAP2_ALIGN(
        ch_ubam_for_alignment,
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
    ubam     = ch_cat_ubam
    bam      = MINIMAP2_ALIGN.out.bam
    index    = MINIMAP2_ALIGN.out.index
    stats    = BAM_STATS_SAMTOOLS.out.stats
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats
}
