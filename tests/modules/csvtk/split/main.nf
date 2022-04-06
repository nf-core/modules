#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CSVTK_SPLIT } from '../../../../modules/csvtk/split/main.nf'

workflow test_csvtk_split_tsv {

    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['generic']['tsv']['test_tsv'], checkIfExists: true) ]
    ]
    in_format = "tsv"
    out_format = "tsv"
    CSVTK_SPLIT ( input, in_format, out_format )
}

workflow test_csvtk_split_csv {

    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['generic']['csv']['test_csv'], checkIfExists: true) ]
    ]
    in_format = "csv"
    out_format = "csv"
    CSVTK_SPLIT( input, in_format, out_format )
}
