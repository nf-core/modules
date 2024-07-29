#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { YAK_COUNT } from '../../../../../modules/nf-core/yak/count/main.nf'

workflow test_yak_count_se {

    input_fq = [ [ id:'test', single_end:true ], 
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
            ] // meta map
    YAK_COUNT ( input_fq )
}

workflow test_yak_count_pe {

    input_fq = [ [ id:'test', single_end:false ], 
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
            ] // meta map
    YAK_COUNT ( input_fq )
}