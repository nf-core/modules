#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MEGAHIT } from '../../../modules/megahit/main.nf'

workflow test_megahit {

    input = [ 
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ] 
    ]

    MEGAHIT ( input )
}

workflow test_megahit_single {

    input = [ 
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    MEGAHIT ( input )
}
