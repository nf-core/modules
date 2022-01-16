#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NEXTCLADE_DATASETGET } from '../../../../modules/nextclade/datasetget/main.nf'
include { NEXTCLADE_RUN        } from '../../../../modules/nextclade/run/main.nf'

workflow test_nextclade_run {
    
    NEXTCLADE_DATASETGET ( 'sars-cov-2' )

    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    NEXTCLADE_RUN ( input, NEXTCLADE_DATASETGET.out.dataset )
}
