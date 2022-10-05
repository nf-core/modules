#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCFLIB_VCFBREAKMULTI } from '../../../../../modules/nf-core/vcflib/vcfbreakmulti/main.nf'

workflow test_vcflib_vcfbreakmulti {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)
    ]

    VCFLIB_VCFBREAKMULTI ( input )
}
