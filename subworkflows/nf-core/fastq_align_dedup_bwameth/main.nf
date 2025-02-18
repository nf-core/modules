include { BWAMETH_ALIGN                                 } from '../../../modules/nf-core/bwameth/align/main'
include { PARABRICKS_FQ2BAMMETH                         } from '../../../modules/nf-core/parabricks/fq2bammeth/main'
include { SAMTOOLS_SORT                                 } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGNMENTS   } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT                             } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_STATS                                } from '../../../modules/nf-core/samtools/stats/main'
include { PICARD_MARKDUPLICATES                         } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUPLICATED } from '../../../modules/nf-core/samtools/index/main'
include { METHYLDACKEL_EXTRACT                          } from '../../../modules/nf-core/methyldackel/extract/main'
include { METHYLDACKEL_MBIAS                            } from '../../../modules/nf-core/methyldackel/mbias/main'

workflow FASTQ_ALIGN_DEDUP_BWAMETH {

    take:
    ch_reads             // channel: [ val(meta), [ reads ] ]
    ch_fasta             // channel: [ val(meta), [ fasta ] ]
    ch_fasta_index       // channel: [ val(meta), [ fasta index ] ]
    ch_bwameth_index     // channel: [ val(meta), [ bwameth index ] ]
    skip_deduplication   // boolean: whether to deduplicate alignments

    main:

    ch_alignment                     = Channel.empty()
    ch_alignment_index               = Channel.empty()
    ch_samtools_flagstat             = Channel.empty()
    ch_samtools_stats                = Channel.empty()
    ch_methydackel_extract_bedgraph  = Channel.empty()
    ch_methydackel_extract_methylkit = Channel.empty()
    ch_methydackel_mbias             = Channel.empty()
    ch_picard_metrics                = Channel.empty()
    ch_multiqc_files                 = Channel.empty()
    ch_versions                      = Channel.empty()

    /*
     * Align with bwameth
     */
    if (params.use_gpu) {
        /*
        * Align with parabricks GPU enabled fq2bammeth implementation of bwameth
        */
        PARABRICKS_FQ2BAMMETH (
            ch_reads,
            ch_fasta,
            ch_bwameth_index,
            [] // known sites
        )
        ch_alignment = PARABRICKS_FQ2BAMMETH.out.bam
        ch_versions  = ch_versions.mix(PARABRICKS_FQ2BAMMETH.out.versions)
    } else {
        /*
        * Align with CPU version of bwameth
        */
        BWAMETH_ALIGN (
            ch_reads,
            ch_fasta,
            ch_bwameth_index
        )
        ch_alignment = BWAMETH_ALIGN.out.bam
        ch_versions  = BWAMETH_ALIGN.out.versions
    }

    /*
     * Sort raw output BAM
     */
    SAMTOOLS_SORT (
        ch_alignment,
        [[:],[]] // [ [meta], [fasta]]
    )
    ch_alignment = SAMTOOLS_SORT.out.bam
    ch_versions  = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    /*
     * Run samtools index on alignment
     */
    SAMTOOLS_INDEX_ALIGNMENTS (
        ch_alignment
    )
    ch_alignment_index = SAMTOOLS_INDEX_ALIGNMENTS.out.bai
    ch_versions        = ch_versions.mix(SAMTOOLS_INDEX_ALIGNMENTS.out.versions)

    /*
     * Run samtools flagstat
     */
    SAMTOOLS_FLAGSTAT (
        ch_alignment.join(ch_alignment_index)
    )
    ch_samtools_flagstat = SAMTOOLS_FLAGSTAT.out.flagstat
    ch_versions          = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    /*
     * Run samtools stats
     */
    SAMTOOLS_STATS (
        ch_alignment.join(ch_alignment_index),
        [[:],[]] // [ [meta], [fasta]]
    )
    ch_samtools_stats = SAMTOOLS_STATS.out.stats
    ch_versions       = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    if (!skip_deduplication) {
        /*
        * Run Picard MarkDuplicates
        */
        PICARD_MARKDUPLICATES (
            ch_alignment,
            ch_fasta,
            ch_fasta_index
        )
        /*
         * Run samtools index on deduplicated alignment
        */
        SAMTOOLS_INDEX_DEDUPLICATED (
            PICARD_MARKDUPLICATES.out.bam
        )
        ch_alignment       = PICARD_MARKDUPLICATES.out.bam
        ch_alignment_index = SAMTOOLS_INDEX_DEDUPLICATED.out.bai
        ch_picard_metrics  = PICARD_MARKDUPLICATES.out.metrics
        ch_versions        = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
        ch_versions        = ch_versions.mix(SAMTOOLS_INDEX_DEDUPLICATED.out.versions)
    }

    /*
     * Extract per-base methylation and plot methylation bias
     */

    METHYLDACKEL_EXTRACT (
        ch_alignment.join(ch_alignment_index),
        ch_fasta.map{ meta, fasta_file -> fasta_file },
        ch_fasta_index.map{ meta, fasta_index -> fasta_index }
    )
    ch_methydackel_extract_bedgraph  = METHYLDACKEL_EXTRACT.out.bedgraph
    ch_methydackel_extract_methylkit = METHYLDACKEL_EXTRACT.out.methylkit
    ch_versions                      = ch_versions.mix(METHYLDACKEL_EXTRACT.out.versions)

    METHYLDACKEL_MBIAS (
        ch_alignment.join(ch_alignment_index),
        ch_fasta.map{ meta, fasta_file -> fasta_file },
        ch_fasta_index.map{ meta, fasta_index -> fasta_index }
    )
    ch_methydackel_mbias = METHYLDACKEL_MBIAS.out.txt
    ch_versions          = ch_versions.mix(METHYLDACKEL_MBIAS.out.versions)

    /*
     * Collect MultiQC inputs
     */
    ch_multiqc_files = ch_picard_metrics.collect{ meta, metrics -> metrics }
                        .mix(ch_samtools_flagstat.collect{ meta, flagstat -> flagstat })
                        .mix(ch_samtools_stats.collect{ meta, stats -> stats  })
                        .mix(ch_methydackel_extract_bedgraph.collect{ meta, bedgraph -> bedgraph  })
                        .mix(ch_methydackel_mbias.collect{ meta, txt -> txt  })

    emit:
    bam                           = ch_alignment                     // channel: [ val(meta), [ bam ]       ]
    bai                           = ch_alignment_index               // channel: [ val(meta), [ bai ]       ]
    samtools_flagstat             = ch_samtools_flagstat             // channel: [ val(meta), [ flagstat ]  ]
    samtools_stats                = ch_samtools_stats                // channel: [ val(meta), [ stats ]     ]
    methydackel_extract_bedgraph  = ch_methydackel_extract_bedgraph  // channel: [ val(meta), [ bedgraph ]  ]
    methydackel_extract_methylkit = ch_methydackel_extract_methylkit // channel: [ val(meta), [ methylkit ] ]
    methydackel_mbias             = ch_methydackel_mbias             // channel: [ val(meta), [ mbias ]     ]
    picard_metrics                = ch_picard_metrics                // channel: [ val(meta), [ metrics ]   ]
    multiqc                       = ch_multiqc_files                 // channel: [ *{html,txt}              ]
    versions                      = ch_versions                      // channel: [ versions.yml             ]
}
