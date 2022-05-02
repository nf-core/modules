#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO } from '../../../modules/busco/main.nf'

workflow test_busco {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
    ]

    BUSCO ( input )
}
