#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCFANNO } from '../../../modules/vcfanno/main.nf'

workflow test_vcfanno {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)
    ]

    toml = file("https://raw.githubusercontent.com/nf-core/test-datasets/8fbd9f99a2feb3f9e39cd3bcdc4a9176a5835673/data/delete_me/vcfanno.toml", 
                checkIfExists: true)

    VCFANNO ( input, toml )
}
