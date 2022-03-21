#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_MARKDUP } from '../../../../modules/samtools/markdup/main.nf'

workflow test_samtools_markdup {
    
    input = [ 
        [ id:'dImpGla2', single_end:false ], // meta map
        file("/lustre/scratch123/tol/teams/tolit/users/ps22/modules/test_data/positionsort.bam", checkIfExists: true) 
    ]

    SAMTOOLS_MARKDUP ( input )
}
