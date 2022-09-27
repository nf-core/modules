#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCFLIB_VCFUNIQ } from '../../../../modules/vcflib/vcfuniq/main.nf'

workflow test_vcflib_vcfuniq {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)
    ]

    VCFLIB_VCFUNIQ ( input )
}
