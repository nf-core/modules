#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAMTOOLS_CONVERT } from '../../../../modules/bamtools/convert/main.nf'

workflow test_bamtools_convert {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
    ]

    BAMTOOLS_CONVERT ( input )
}
