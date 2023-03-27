include { HISAT2_ALIGN            } from '../../../modules/nf-core/hisat2/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_HISAT2 {

    take:
    ch_reads       // channel: [ val(meta), path(reads) ]
    ch_index       // channel: [ path(index) ]      (/path/to/hisat2/index)
    ch_splicesites // channel: [ path(txt) ]        (/path/to/genome.splicesites.txt)
    ch_fasta    // channel: [ path(fasta) ]

    main:

    ch_versions = Channel.empty()


    //
    // Map reads with HISAT2
    //
    HISAT2_ALIGN ( ch_reads, ch_index, ch_splicesites )
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( HISAT2_ALIGN.out.bam, ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)


    emit:
    orig_bam = HISAT2_ALIGN.out.bam                 // channel: [ val(meta), path(bam) ]
    summary  = HISAT2_ALIGN.out.summary             // channel: [ val(meta), path(log) ]
    fastq    = HISAT2_ALIGN.out.fastq               // channel: [ val(meta), path(fastq)]

    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), path(bam) ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), path(bai) ]
    csi      = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), path(csi) ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                          // channel: [ path(versions.yml) ]
}

