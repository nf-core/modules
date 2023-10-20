//
// Alignment with Bowtie2
//

include { BOWTIE2_ALIGN           } from '../../../modules/nf-core/bowtie2/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_BOWTIE2 {
    take:
    ch_reads            // channel: [ val(meta), path(reads) ]
    ch_index            // channel: [ path(index) ]     (/path/to/bowtie2/index/)
    ch_save_unaligned   // channel: [ val(unaligned)]
    ch_sort_bam         // channel: [ val(bam)]
    ch_fasta            // channel:[ path(fasta) ]     (/path/to/reference.fasta)

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with Bowtie2
    //
    BOWTIE2_ALIGN ( ch_reads, ch_index, ch_save_unaligned, ch_sort_bam )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( BOWTIE2_ALIGN.out.bam, ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam_orig         = BOWTIE2_ALIGN.out.bam          // channel: [ val(meta), val(bam) ]
    log_out          = BOWTIE2_ALIGN.out.log          // channel: [ val(meta), path(log) ]
    fastq            = BOWTIE2_ALIGN.out.fastq        // channel: [ val(meta), path(fastq) ]

    bam              = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), val(bam) ]
    bai              = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), path(bai) ]
    csi              = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), path(csi) ]
    stats            = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat         = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats         = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]

    versions         = ch_versions                      // channel: [ path(versions.yml) ]
}
