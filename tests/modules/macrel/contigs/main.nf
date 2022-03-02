#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MACREL_CONTIGS } from '../../../../modules/macrel/contigs/main.nf'

workflow test_macrel_contigs {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true) 
    ]

    MACREL_CONTIGS ( input )
}
