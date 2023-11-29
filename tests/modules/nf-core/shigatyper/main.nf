#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHIGATYPER } from '../../../../modules/nf-core/shigatyper/main.nf'

workflow test_shigatyper_pe {
    
    input = [ 
        [ id:'test', single_end:false, is_ont:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
          file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
    ]

    SHIGATYPER ( input )
}

workflow test_shigatyper_se {
    
    input = [ 
        [ id:'test', single_end:true, is_ont:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
    ]

    SHIGATYPER ( input )
}

workflow test_shigatyper_ont {
    
    input = [ 
        [ id:'test', single_end:true, is_ont:true ], // meta map
        [ file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ]
    ]

    SHIGATYPER ( input )
}
