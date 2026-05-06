/*
 * Map reads, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { CHROMAP_CHROMAP         } from '../../../modules/nf-core/chromap/chromap/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_CHROMAP {
    take:
    ch_reads // channel (mandatory): [ val(meta), [ reads ] ]
    ch_index // channel (mandatory): [ val(meta2, [ index ] ]
    ch_fasta_fai // channel (mandatory): [ val(meta2, [ fasta ], [ fai ] ]
    ch_barcodes // channel (optional):  [ barcodes ]
    ch_whitelist // channel (optional):  [ whitelist ]
    ch_chr_order // channel (optional):  [ chr_order ]
    ch_pairs_chr_order // channel (optional):  [ pairs_chr_order ]

    main:

    //
    // Remap ch_fasta_fai to ch_fasta
    //
    ch_fasta = ch_fasta_fai.map { meta, fasta , _fai -> [ meta, fasta ] }

    //
    // Map reads with CHROMAP
    //
    CHROMAP_CHROMAP(ch_reads, ch_fasta, ch_index, ch_barcodes, ch_whitelist, ch_chr_order, ch_pairs_chr_order)

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS(CHROMAP_CHROMAP.out.bam, ch_fasta_fai)

    emit:
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam // channel: [ val(meta), [ bam ] ]
    index    = BAM_SORT_STATS_SAMTOOLS.out.index // channel: [ val(meta), [ index ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
}
