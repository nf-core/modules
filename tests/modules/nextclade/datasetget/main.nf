#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NEXTCLADE_DATASETGET } from '../../../../modules/nextclade/datasetget/main.nf'

workflow test_nextclade_datasetget {

    NEXTCLADE_DATASETGET ( 'sars-cov-2' )
}
