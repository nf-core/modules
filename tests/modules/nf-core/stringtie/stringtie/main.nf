#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STRINGTIE_STRINGTIE } from '../../../../../modules/nf-core/stringtie/stringtie/main.nf'

//
// Test with forward strandedness
//
workflow test_stringtie_forward {
    input = [
        [ id:'test', strandedness:'forward' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
    ]

    STRINGTIE_STRINGTIE ( input, [] )
}

//
// Test with forward strandedness and reference annotation
//
workflow test_stringtie_forward_annotation {
    input = [
        [ id:'test', strandedness:'forward' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
    ]
    annotation_gtf = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)

    STRINGTIE_STRINGTIE ( input, annotation_gtf )
}

//
// Test with reverse strandedness
//
workflow test_stringtie_reverse {
    input = [
        [ id:'test', strandedness:'reverse' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
    ]

    STRINGTIE_STRINGTIE ( input, [] )
}

//
// Test with reverse strandedness and reference annotation
//
workflow test_stringtie_reverse_annotation {
    input = [
        [ id:'test', strandedness:'reverse' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
    ]
    annotation_gtf = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)

    STRINGTIE_STRINGTIE ( input, annotation_gtf )
}
