#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_BBNORM } from '../../../../../modules/nf-core/bbmap/bbnorm/main.nf'

workflow test_bbmap_bbnorm_se {
    
    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    BBMAP_BBNORM ( input )
}

workflow test_bbmap_bbnorm_pe {

    input = [ 
        [ id:'test', single_end:false ], // meta map
        [  
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) 
        ]
    ]

    BBMAP_BBNORM ( input )
}

workflow test_bbmap_bbnorm_interleaved {

    input  = [
        [id:'test', single_end:true ],
        [ file(params.test_data['sarscov2']['illumina']['test_interleaved_fastq_gz'], checkIfExists: true)]
    ]

    BBMAP_BBNORM ( input )
}

 workflow test_bbmap_bbnorm_multiple_input {
     input = [
        [id:'test', single_end:true ],
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    BBMAP_BBNORM ( input )
}
