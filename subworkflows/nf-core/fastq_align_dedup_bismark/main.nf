include { BISMARK_ALIGN                                 } from '../../../modules/nf-core/bismark/align/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_DEDUPLICATED   } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUPLICATED } from '../../../modules/nf-core/samtools/index/main'
include { BISMARK_DEDUPLICATE                           } from '../../../modules/nf-core/bismark/deduplicate/main'
include { BISMARK_METHYLATIONEXTRACTOR                  } from '../../../modules/nf-core/bismark/methylationextractor/main'
include { BISMARK_COVERAGE2CYTOSINE                     } from '../../../modules/nf-core/bismark/coverage2cytosine/main'
include { BISMARK_REPORT                                } from '../../../modules/nf-core/bismark/report/main'
include { BISMARK_SUMMARY                               } from '../../../modules/nf-core/bismark/summary/main'

workflow FASTQ_ALIGN_DEDUP_BISMARK {

    take:
    ch_reads             // channel: [ val(meta), [ reads ] ]
    ch_fasta             // channel: [ val(meta), [ fasta ] ]
    ch_bismark_index     // channel: [ val(meta), [ bismark index ] ]
    skip_deduplication   // boolean: whether to deduplicate alignments
    cytosine_report      // boolean: whether the run coverage2cytosine

    main:

    ch_versions = Channel.empty()

    /*
     * Align with bismark
     */
    BISMARK_ALIGN (
        ch_reads,
        ch_fasta,
        ch_bismark_index
    )

    ch_versions = ch_versions.mix(BISMARK_ALIGN.out.versions)


    if (skip_deduplication) {
        alignments        = BISMARK_ALIGN.out.bam
        alignment_reports = BISMARK_ALIGN.out.report.map{ meta, report -> [ meta, report, [] ] }
    } else {
        /*
        * Run deduplicate_bismark
        */
        BISMARK_DEDUPLICATE(BISMARK_ALIGN.out.bam)

        alignments        = BISMARK_DEDUPLICATE.out.bam
        alignment_reports = BISMARK_ALIGN.out.report.join(BISMARK_DEDUPLICATE.out.report)
        ch_versions       = ch_versions.mix(BISMARK_DEDUPLICATE.out.versions)
    }

    /*
     * Run bismark_methylation_extractor
     */
    BISMARK_METHYLATIONEXTRACTOR (
        alignments,
        ch_bismark_index
    )

    ch_versions = ch_versions.mix(BISMARK_METHYLATIONEXTRACTOR.out.versions)

    /*
     * Run coverage2cytosine
     */
    if (cytosine_report) {
        BISMARK_COVERAGE2CYTOSINE (
            BISMARK_METHYLATIONEXTRACTOR.out.coverage,
            ch_fasta,
            ch_bismark_index
        )
        ch_cytosine_report_summary = BISMARK_COVERAGE2CYTOSINE.out.report.collect{ it[1] }
                                        .mix(BISMARK_COVERAGE2CYTOSINE.out.summary.collect{ it[1] }.ifEmpty([]))
        ch_versions                = ch_versions.mix(BISMARK_COVERAGE2CYTOSINE.out.versions)
    }

    /*
     * Generate bismark sample reports
     */
    BISMARK_REPORT (
        alignment_reports
            .join(BISMARK_METHYLATIONEXTRACTOR.out.report)
            .join(BISMARK_METHYLATIONEXTRACTOR.out.mbias)
    )

    ch_versions = ch_versions.mix(BISMARK_REPORT.out.versions)

    /*
     * Generate bismark summary report
     */
    BISMARK_SUMMARY (
        BISMARK_ALIGN.out.bam.collect{ it[1].name }.ifEmpty([]),
        alignment_reports.collect{ it[1] }.ifEmpty([]),
        alignment_reports.collect{ it[2] }.ifEmpty([]),
        BISMARK_METHYLATIONEXTRACTOR.out.report.collect{ it[1] }.ifEmpty([]),
        BISMARK_METHYLATIONEXTRACTOR.out.mbias.collect{ it[1] }.ifEmpty([])
    )

    ch_versions = ch_versions.mix(BISMARK_SUMMARY.out.versions)

    /*
     * MODULE: Run samtools sort on dedup bam
     */
    SAMTOOLS_SORT_DEDUPLICATED (
        alignments,
        [[:],[]] // Empty map and list as is optional input but required for nextflow
    )

    ch_versions = ch_versions.mix(SAMTOOLS_SORT_DEDUPLICATED.out.versions)

    /*
     * MODULE: Run samtools index on dedup bam
     */
    SAMTOOLS_INDEX_DEDUPLICATED(SAMTOOLS_SORT_DEDUPLICATED.out.bam)

    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_DEDUPLICATED.out.versions)

    /*
     * Collect MultiQC inputs
     */
    BISMARK_SUMMARY.out.summary.ifEmpty([])
        .mix(alignment_reports.collect{ it[1] })
        .mix(alignment_reports.collect{ it[2] })
        .mix(BISMARK_METHYLATIONEXTRACTOR.out.report.collect{ it[1] })
        .mix(BISMARK_METHYLATIONEXTRACTOR.out.mbias.collect{ it[1] })
        .mix(BISMARK_REPORT.out.report.collect{ it[1] })
        .set{ multiqc_files }

    emit:
    dedup           = SAMTOOLS_SORT_DEDUPLICATED.out.bam // channel: [ val(meta), [ bam ] ] ## sorted, possibly deduplicated BAM
    cytosine_report = ch_cytosine_report_summary         // channel: [ val(meta), [ report summary ] ] ## bismark cytosine report and summary
    multiqc         = multiqc_files                      // path: *{html,txt}
    versions        = ch_versions                        // path: *.version.txt
}

