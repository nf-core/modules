//
// Alignment with Bowtie2
//

params.align_options          = [:]
params.samtools_sort_options  = [:]
params.samtools_index_options = [:]
params.samtools_stats_options = [:]

include { BOWTIE2_ALIGN     } from '../../../modules/bowtie2/align/main' addParams( options: params.align_options    )
include { BAM_SORT_SAMTOOLS } from '../bam_sort_samtools/main'           addParams( sort_options: params.samtools_sort_options, index_options: params.samtools_index_options, stats_options: params.samtools_stats_options )

workflow ALIGN_BOWTIE2 {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index // channel: /path/to/bowtie2/index/

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with Bowtie2
    //
    BOWTIE2_ALIGN ( reads, index )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_SAMTOOLS ( BOWTIE2_ALIGN.out.bam )
    ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)

    emit:
    bam_orig          = BOWTIE2_ALIGN.out.bam          // channel: [ val(meta), bam   ]
    log_out           = BOWTIE2_ALIGN.out.log          // channel: [ val(meta), log   ]
    fastq             = BOWTIE2_ALIGN.out.fastq        // channel: [ val(meta), fastq ]

    bam              = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai              = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    csi              = BAM_SORT_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    stats            = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions       = ch_versions                      // channel: [ versions.yml ]
}
