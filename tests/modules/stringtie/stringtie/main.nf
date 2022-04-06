#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STRINGTIE } from '../../../../modules/stringtie/stringtie/main.nf'
//
// Test with forward strandedness
//
workflow test_stringtie_forward {
    input = [
        [ id:'test', strandedness:'forward' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
    ]
    annotation_gtf = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)

    STRINGTIE ( input, annotation_gtf )
}

//
// Test with reverse strandedness
//
workflow test_stringtie_reverse {
    input = [
        [ id:'test', strandedness:'reverse' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
    ]
    annotation_gtf = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)

    STRINGTIE ( input, annotation_gtf )
}
