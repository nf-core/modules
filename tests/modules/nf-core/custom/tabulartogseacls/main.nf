#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUSTOM_TABULARTOGSEACLS } from '../../../../../modules/nf-core/custom/tabulartogseacls/main.nf'

process csv_to_tsv {

    input:
    path expression_matrix

    output:
    path 'test.tsv'

    script:
    """
    sed 's/,/\\t/g' ${expression_matrix} > test.tsv.tmp
    mv test.tsv.tmp test.tsv
    """

}
workflow test_custom_tabulartogseacls {
    
    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    
    input = [
        [ id:'test', variable:'treatment' ], // meta map
        expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    ]

    CUSTOM_TABULARTOGSEACLS ( input )
}

workflow test_custom_tabulartogseacls_tsv {
    
    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    
    input = csv_to_tsv(file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true))
        .map{
            tuple([ id:'test', variable:'treatment' ], it)
        }

    CUSTOM_TABULARTOGSEACLS ( input )
}
