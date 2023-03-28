include { IVAR_TRIM               } from '../../../modules/nf-core/ivar/trim/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow BAM_TRIM_PRIMERS_IVAR {
    take:
    ch_bam   // channel: [ val(meta), path(bam), path(bai) ]
    ch_bed   // channel: [ path(bed) ]
    ch_fasta // channel: [ path(fasta) ]

    main:
    ch_versions = Channel.empty()

    //
    // iVar trim primers
    //
    IVAR_TRIM (
        ch_bam,
        ch_bed
    )
    ch_versions = ch_versions.mix(IVAR_TRIM.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS (
        IVAR_TRIM.out.bam,
        ch_fasta
    )

    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions.first())

    emit:
    bam_orig = IVAR_TRIM.out.bam                    // channel: [ val(meta), path(bam) ]
    log_out  = IVAR_TRIM.out.log                    // channel: [ val(meta), path(log) ]

    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), path(bam) ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), path(bai) ]
    csi      = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), path(csi) ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                          // channel: [ path(versions.yml) ]
}
