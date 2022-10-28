/*
 * Map reads, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { CHROMAP_CHROMAP         } from '../../../modules/nf-core/chromap/chromap/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_CHROMAP {
    take:
    ch_reads           // channel: [ val(meta), [ reads ] ]
    ch_index           // channel: [ val(meta2, [ index ] ]
    ch_fasta           // channel: [ val(meta2, [ fasta ] ]
    ch_barcodes        // channel: [ barcodes ]
    ch_whitelist       // channel: [ whitelist ]
    ch_chr_order       // channel: [ chr_order ]
    ch_pairs_chr_order // channel: [ pairs_chr_order ]

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
    BAM_SORT_STATS_SAMTOOLS(CHROMAP_CHROMAP.out.bam, [])
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam               = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai               = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats             = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat          = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats          = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions          = ch_versions                          //    path: versions.yml
}
