//
// Alignment with Bowtie2
//

params.align_options    = [:]
params.samtools_options = [:]

include { BOWTIE2_ALIGN     } from '../../modules/nf-core/software/bowtie2/align/main' addParams( options: params.align_options    )
include { BAM_SORT_SAMTOOLS } from './bam_sort_samtools'                      addParams( options: params.samtools_options )

workflow ALIGN_BOWTIE2 {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index // channel: /path/to/bowtie2/index/

    main:

    //
    // Map reads with Bowtie2
    //
    BOWTIE2_ALIGN ( reads, index )

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_SAMTOOLS ( BOWTIE2_ALIGN.out.bam )

    emit:
    bam_orig         = BOWTIE2_ALIGN.out.bam          // channel: [ val(meta), bam   ]
    log_out          = BOWTIE2_ALIGN.out.log          // channel: [ val(meta), log   ]
    fastq            = BOWTIE2_ALIGN.out.fastq        // channel: [ val(meta), fastq ]
    bowtie2_version  = BOWTIE2_ALIGN.out.version      //    path: *.version.txt

    bam              = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai              = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats            = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    samtools_version = BAM_SORT_SAMTOOLS.out.version  //    path: *.version.txt
}
