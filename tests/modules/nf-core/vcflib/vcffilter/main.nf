#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCFLIB_VCFFILTER } from '../../../../../modules/nf-core/vcflib/vcffilter/main.nf'

workflow test_vcflib_vcffilter {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)
    ]

    VCFLIB_VCFFILTER ( input )
}

