#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GRAPHTYPER_VCFCONCATENATE } from '../../../../../modules/nf-core/graphtyper/vcfconcatenate/main.nf'

workflow test_graphtyper_vcfconcatenate {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true)
        ]
    ]

    GRAPHTYPER_VCFCONCATENATE ( input )
}
