//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
//

include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'

workflow BAM_SORT_STATS_SAMTOOLS {
    take:
    ch_bam   // channel: [ val(meta), [ bam ] ]
    ch_fasta // channel: [ val(meta), path(fasta) ]

    main:
    SAMTOOLS_SORT ( ch_bam, ch_fasta, 'bai' )

    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_SORT.out.bai, by: [0], remainder: true)
        .set { ch_bam_bai }

    BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta )

    emit:
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_SORT.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_SORT.out.csi          // channel: [ val(meta), [ csi ] ]

    stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
}
