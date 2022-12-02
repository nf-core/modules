#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUSTOM_TABULARTOGSEAGCT } from '../../../../../modules/nf-core/custom/tabulartogseagct/main.nf'

process tsv_to_csv {

    input:
    path expression_matrix

    output:
    path 'test.csv'

    script:
    """
    sed 's/\t/,/g' ${expression_matrix} > test.csv.tmp
    mv test.csv.tmp test.csv
    """

}

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

workflow test_custom_tabulartogseagct {
    
    input = strip_column(file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true))
        .map{
            tuple([ id:'test' ], it)
        }

    CUSTOM_TABULARTOGSEAGCT ( input )
}

workflow test_custom_tabulartogseagct_csv {
    
    strip_column(file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true))
    
    input = tsv_to_csv(strip_column.out)    
        .map{
            tuple([ id:'test' ], it)
        }

    CUSTOM_TABULARTOGSEAGCT ( input )
}
