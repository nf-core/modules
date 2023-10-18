#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CLAME } from '../../../../modules/nf-core/clame/main.nf'
include { GUNZIP } from '../../../../modules/nf-core/gunzip/main.nf'

workflow test_clame {
    
   input = [ 
        [ id:'test1', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true)
    ]
    
    GUNZIP ( input )

    CLAME ( GUNZIP.out.gunzip )
}
