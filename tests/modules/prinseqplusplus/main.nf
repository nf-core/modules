#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRINSEQPLUSPLUS } from '../../../modules/prinseqplusplus/main.nf'

workflow test_prinseqplusplus_single_end {

        input = [ [ id:'test', single_end:true ], // meta map
            [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
        ]

    PRINSEQPLUSPLUS ( input )
}

workflow test_prinseqplusplus_paired_end {

    input = [ [ id:'test', single_end:false ], // meta map
            [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
        ]

    PRINSEQPLUSPLUS ( input )
}
