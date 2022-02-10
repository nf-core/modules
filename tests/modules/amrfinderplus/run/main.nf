#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AMRFINDERPLUS_UPDATE } from '../../../../modules/amrfinderplus/update/main.nf'
include { AMRFINDERPLUS_RUN } from '../../../../modules/amrfinderplus/run/main.nf'

workflow test_amrfinderplus_run {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['haemophilus_influenzae']['genome']['genome_fna_gz'], checkIfExists: true) 
    ]

    AMRFINDERPLUS_UPDATE ( )
    AMRFINDERPLUS_RUN ( input, AMRFINDERPLUS_UPDATE.out.db )
}
