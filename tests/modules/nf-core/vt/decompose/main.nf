#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VT_DECOMPOSE } from '../../../../../modules/nf-core/vt/decompose/main.nf'

workflow test_vt_decompose {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        []
    ]

    VT_DECOMPOSE ( input )
}

workflow test_vt_decompose_intervals {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    ]

    VT_DECOMPOSE ( input )
}
