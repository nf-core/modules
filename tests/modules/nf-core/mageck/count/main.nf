#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAGECK_COUNT } from '../../../../../modules/nf-core/mageck/count/main.nf'

workflow test_mageck_count_fastq {
    def input = []
    input = [ [ id:'test,test2', single_end:true] , // meta map
        [ file(params.test_data['mus_musculus']['illumina']['test_1_fastq_gz'], checkIfExists: true),
          file(params.test_data['mus_musculus']['illumina']['test_2_fastq_gz'], checkIfExists: true)]
    ]
    library = file(params.test_data['mus_musculus']['csv']['library'], checkIfExists: true)

    MAGECK_COUNT ( input, library )
}


workflow test_mageck_count_counts {
    def input = []
    input =  [[ id:'test', single_end:true] , // meta map
        [file(params.test_data['mus_musculus']['csv']['count_table'], checkIfExists: true)]]
    library = file(params.test_data['mus_musculus']['csv']['library'], checkIfExists: true)

    MAGECK_COUNT ( input, library)
}
