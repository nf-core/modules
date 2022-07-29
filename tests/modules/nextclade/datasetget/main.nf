#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NEXTCLADE_DATASETGET } from '../../../../modules/nextclade/datasetget/main.nf'

workflow test_nextclade_datasetget {

    dataset = 'sars-cov-2'
    reference = 'MN908947'
    tag = '2022-01-18T12:00:00Z'

    NEXTCLADE_DATASETGET ( dataset, reference, tag )

}
