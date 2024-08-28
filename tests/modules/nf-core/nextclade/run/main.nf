#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NEXTCLADE_DATASETGET } from '../../../../../modules/nf-core/nextclade/datasetget/main.nf'
include { NEXTCLADE_RUN        } from '../../../../../modules/nf-core/nextclade/run/main.nf'

workflow test_nextclade_run {

    dataset = 'sars-cov-2'
    tag = '2024-01-16--20-31-02Z'

    NEXTCLADE_DATASETGET ( dataset, tag )

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    NEXTCLADE_RUN ( input, NEXTCLADE_DATASETGET.out.dataset )
}

