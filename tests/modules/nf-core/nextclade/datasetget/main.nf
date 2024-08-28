#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NEXTCLADE_DATASETGET } from '../../../../../modules/nf-core/nextclade/datasetget/main.nf'

workflow test_nextclade_datasetget {

    dataset = 'nextstrain/sars-cov-2/wuhan-hu-1/orfs'
    tag = '2024-01-16--20-31-02Z'

    NEXTCLADE_DATASETGET ( dataset, tag )

}

