#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCFLIB_VCFFILTER } from '../../../../../modules/nf-core/vcflib/vcffilter/main.nf'

workflow test_vcflib_vcffilter_vcf {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    VCFLIB_VCFFILTER ( input, [[],[]] )
}

workflow test_vcflib_vcffilter_gz {

    vcf_input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true)
    ]
    tbi_input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)
    ]

    VCFLIB_VCFFILTER ( vcf_input, tbi_input )
}

