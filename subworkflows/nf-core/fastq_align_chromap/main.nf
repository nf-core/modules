/*
 * Map reads, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { CHROMAP_CHROMAP         } from '../../../modules/nf-core/chromap/chromap/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_CHROMAP {
    take:
    ch_reads           // channel (mandatory): [ val(meta), path(reads) ]
    ch_index           // channel (mandatory): [ val(meta), path(index) ]
    ch_fasta           // channel (mandatory): [ val(meta), path(fasta) ]
    ch_barcodes        // channel (optional):  [ path(barcodes) ]
    ch_whitelist       // channel (optional):  [ path(whitelist) ]
    ch_chr_order       // channel (optional):  [ path(chr_order) ]
    ch_pairs_chr_order // channel (optional):  [ path(pairs_chr_order) ]

    main:
    ch_versions = Channel.empty()

    //
    // Map reads with CHROMAP
    //
    CHROMAP_CHROMAP(ch_reads, ch_fasta, ch_index, ch_barcodes, ch_whitelist, ch_chr_order, ch_pairs_chr_order)
    ch_versions = ch_versions.mix(CHROMAP_CHROMAP.out.versions)

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS(CHROMAP_CHROMAP.out.bam, ch_fasta.map { it[1] })
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam               = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), path(bam) ]
    bai               = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), path(bai) ]
    stats             = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat          = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats          = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]

    versions          = ch_versions                          // channel: [(path: versions.yml)]
}
3
