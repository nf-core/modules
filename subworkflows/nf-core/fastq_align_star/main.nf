include { STAR_ALIGN              } from '../../../modules/nf-core/star/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_STAR {

    take:
    ch_reads                    // channel: [ val(meta), [ path(reads) ] ]
    ch_index                    // channel: [ path(index) ]
    ch_gtf                      // channel: [ path(gtf) ]
    val_star_ignore_sjdbgtf     // boolean: when using pre-built STAR indices do not re-extract and use splice junctions from the GTF file
    val_seq_platform            // string : sequencing platform
    val_seq_center              // string : sequencing center
    ch_fasta                    // channel: [ val(meta), path(fasta) ]

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with STAR
    //
    STAR_ALIGN ( ch_reads, ch_index, ch_gtf, val_star_ignore_sjdbgtf, val_seq_platform, val_seq_center )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( STAR_ALIGN.out.bam, ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:

    orig_bam       = STAR_ALIGN.out.bam             // channel: [ val(meta), path(bam)            ]
    log_final      = STAR_ALIGN.out.log_final       // channel: [ val(meta), path(log_final)      ]
    log_out        = STAR_ALIGN.out.log_out         // channel: [ val(meta), path(log_out)        ]
    log_progress   = STAR_ALIGN.out.log_progress    // channel: [ val(meta), path(log_progress)   ]
    bam_sorted     = STAR_ALIGN.out.bam_sorted      // channel: [ val(meta), path(bam)            ]
    bam_transcript = STAR_ALIGN.out.bam_transcript  // channel: [ val(meta), path(bam)            ]
    fastq          = STAR_ALIGN.out.fastq           // channel: [ val(meta), path(fastq)          ]
    tab            = STAR_ALIGN.out.tab             // channel: [ val(meta), path(tab)            ]

    bam            = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), path(bam) ]
    bai            = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), path(bai) ]
    stats          = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat       = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats       = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]

    versions       = ch_versions                    // channel: [ path(versions.yml) ]
}
