include { FASTQ_ALIGN_BWA                                   } from '../fastq_align_bwa/main'
include { PICARD_ADDORREPLACEREADGROUPS                     } from '../../../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_MARKDUPLICATES                             } from '../../../modules/nf-core/picard/markduplicates/main'  
include { SAMTOOLS_INDEX                                    } from '../../../modules/nf-core/samtools/index/main'

workflow FASTQ_ALIGN_DEDUP_BWAMEM {

    take:
    ch_reads             // channel: [ val(meta), [ reads ] ]
    ch_fasta             // channel: [ val(meta), [ fasta ] ]
    ch_fasta_index       // channel: [ val(meta), [ fasta index ] ]
    ch_bwamem_index      // channel: [ val(meta), [ bwam index ] ]
    skip_deduplication   // boolean: whether to deduplicate alignments
    use_gpu              // boolean: whether to use GPU or CPU for bwamem alignment

    main:

    ch_alignment                     = Channel.empty()
    ch_alignment_index               = Channel.empty()
    ch_flagstat                      = Channel.empty()
    ch_stats                         = Channel.empty()
    ch_picard_metrics                = Channel.empty()
    ch_multiqc_files                 = Channel.empty()
    ch_versions                      = Channel.empty()

    FASTQ_ALIGN_BWA (
        ch_reads,
        ch_bwamem_index,
        true, // val_sort_bam hardcoded to true
        ch_fasta
    )
    ch_alignment        = ch_alignment.mix(FASTQ_ALIGN_BWA.out.bam)
    ch_alignment_index  = ch_alignment.mix(FASTQ_ALIGN_BWA.out.bai)
    ch_stats            = ch_alignment.mix(FASTQ_ALIGN_BWA.out.stats)    // channel: [ val(meta), path(stats) ]
    ch_flagstat         = ch_alignment.mix(FASTQ_ALIGN_BWA.out.flagstat) // channel: [ val(meta), path(flagstat) ]
    ch_idxstats         = ch_alignment.mix(FASTQ_ALIGN_BWA.out.idxstats) // channel: [ val(meta), path(idxstats) ]
    ch_versions         = ch_versions.mix(FASTQ_ALIGN_BWA.out.versions.first())

    if (!skip_deduplication) {
        /*
         * Run Picard AddOrReplaceReadGroups to add read group (RG) to reads in bam file
         */
        PICARD_ADDORREPLACEREADGROUPS (
            ch_alignment,
            ch_fasta,
            ch_fasta_index
        )
        ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS.out.versions.first())

        /*
         * Run Picard MarkDuplicates with the --REMOVE_DUPLICATES true flag
         */

        PICARD_MARKDUPLICATES (
            PICARD_ADDORREPLACEREADGROUPS.out.bam,
            ch_fasta,
            ch_fasta_index
        )
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

        /*
         * Run samtools index on deduplicated alignment
         */
        SAMTOOLS_INDEX (
            PICARD_MARKDUPLICATES.out.bam
        )
        ch_alignment       = PICARD_MARKDUPLICATES.out.bam
        ch_alignment_index = SAMTOOLS_INDEX.out.bai
        ch_picard_metrics  = PICARD_MARKDUPLICATES.out.metrics
        ch_versions        = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    }

    /*
     * Collect MultiQC inputs
     */
    ch_multiqc_files = ch_picard_metrics.collect{ meta, metrics -> metrics }
                        .mix(ch_flagstat.collect{ meta, flagstat -> flagstat })
                        .mix(ch_stats.collect{ meta, stats -> stats  })
                        .mix(ch_idxstats.collect{ meta, stats -> stats  })

    emit:
    bam                           = ch_alignment                     // channel: [ val(meta), [ bam ]       ]
    bai                           = ch_alignment_index               // channel: [ val(meta), [ bai ]       ]
    samtools_flagstat             = ch_flagstat                      // channel: [ val(meta), [ flagstat ]  ]
    samtools_stats                = ch_stats                         // channel: [ val(meta), [ stats ]     ]
    samtools_index_stats          = ch_idxstats                      // channel: [ val(meta), [ idxstats ]  ]
    picard_metrics                = ch_picard_metrics                // channel: [ val(meta), [ metrics ]   ]
    multiqc                       = ch_multiqc_files                 // channel: [ *{html,txt}              ]
    versions                      = ch_versions                      // channel: [ versions.yml             ]
}