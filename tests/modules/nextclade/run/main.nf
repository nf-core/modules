#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NEXTCLADE_DATASETGET } from '../../../../modules/nextclade/datasetget/main.nf'
include { NEXTCLADE_RUN        } from '../../../../modules/nextclade/run/main.nf'

workflow test_nextclade_run {

    dataset = 'sars-cov-2'
    reference = 'MN908947'
    tag = '2022-01-18T12:00:00Z'

    NEXTCLADE_DATASETGET ( dataset, reference, tag )

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    NEXTCLADE_RUN ( input, NEXTCLADE_DATASETGET.out.dataset )
}

