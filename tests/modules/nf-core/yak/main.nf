#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { YAK } from '../../../../modules/nf-core/yak/main.nf'

workflow test_yak {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    YAK ( input )
}
