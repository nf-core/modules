#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STRINGTIE       } from '../../../../modules/stringtie/stringtie/main.nf'
include { STRINGTIE_MERGE } from '../../../../modules/stringtie/merge/main.nf'

/*
 * Test with forward strandedness
 */
workflow test_stringtie_forward_merge {
    input = [
        [ id:'test', strandedness:'forward' ], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
    ]
    annotation_gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    STRINGTIE ( input, annotation_gtf )
    STRINGTIE
        .out
        .transcript_gtf
        .map { it -> it[1] }
        .set { stringtie_gtf }

    STRINGTIE_MERGE ( stringtie_gtf, annotation_gtf )
}

/*
 * Test with reverse strandedness
 */
workflow test_stringtie_reverse_merge {
    input = [
        [ id:'test', strandedness:'reverse' ], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
    ]
    annotation_gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    STRINGTIE ( input, annotation_gtf )
    STRINGTIE
        .out
        .transcript_gtf
        .map { it -> it[1] }
        .set { stringtie_gtf }

    STRINGTIE_MERGE ( stringtie_gtf, annotation_gtf )
}
