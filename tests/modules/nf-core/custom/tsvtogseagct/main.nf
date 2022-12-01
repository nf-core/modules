#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUSTOM_TSVTOGSEAGCT } from '../../../../../modules/nf-core/custom/tsvtogseagct/main.nf'

// Example expression matrix comes from rnaseq workflow, which puts symbols in
// the second column, strip columns

process strip_column {

    input:
    path tsv

    output:
    path 'test.tsv'

    script:
    """
    cut -f2 --complement $tsv > test.tsv
    """
}

workflow test_custom_tsvtogseagct {
    
    input = strip_column(file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true))
        .map{
            tuple([ id:'test' ], it)
        }

    CUSTOM_TSVTOGSEAGCT ( input )
}
