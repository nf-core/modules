#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FILTLONG } from '../../../modules/filtlong/main.nf'

workflow test_filtlong {
    
    input = [ [ id:'test', single_end:false ], // meta map
              [],
              [ file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ]
            ]

    FILTLONG ( input )
}

workflow test_filtlong_illumina_se {
    
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ],
              [ file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ]
            ]

    FILTLONG ( input )
}

workflow test_filtlong_illumina_pe {
    
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ],
              [ file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ]
            ]

    FILTLONG ( input )
}
