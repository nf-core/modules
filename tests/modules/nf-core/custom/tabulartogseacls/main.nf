#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUSTOM_TABULARTOGSEACLS } from '../../../../../modules/nf-core/custom/tabulartogseacls/main.nf'

expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)

workflow test_custom_tabulartogseacls {
    
    input = [
        [ id:'test', variable:'treatment' ], // meta map
        expression_sample_sheet
    ]

    CUSTOM_TABULARTOGSEACLS ( input )
}

workflow test_custom_tabulartogseacls_tsv {
   
    input = Channel.fromPath(expression_sample_sheet)
       .splitCsv(header: false)
       .map{
           it.join('\t')
       }
       .collectFile(name: 'test.tsv', newLine: true, sort: false)
       .map{
           [ [ id:'test', variable:'treatment' ], it]
       }

    CUSTOM_TABULARTOGSEACLS ( input )
}

workflow test_custom_tabulartogseacls_tsv_override {
   
    input = Channel.fromPath(expression_sample_sheet)
       .splitCsv(header: false)
       .map{
           it.join('\t')
       }
       .collectFile(name: 'test.csv', newLine: true, sort: false)
       .map{
           [ [ id:'test', variable:'treatment' ], it]
       }

    CUSTOM_TABULARTOGSEACLS ( input )
}
