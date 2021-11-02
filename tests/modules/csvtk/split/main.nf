#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CSVTK_SPLIT as CSVTK_SPLIT_TSV_TSV } from '../../../../modules/csvtk/split/main.nf' addParams( options: [args: "-C '&' --fields 'first_name' "])
include { CSVTK_SPLIT as CSVTK_SPLIT_CSV_CSV } from '../../../../modules/csvtk/split/main.nf' addParams( options: [args: "-C '&' --fields 'first_name' "])

workflow test_csvtk_split_tsv_tsv {

    input = [
        [ id:'test' ], // meta map
        [ file("https://raw.githubusercontent.com/SusiJo/test-datasets/modules/data/generic/tsv/test.tsv", checkIfExists: true) ]
    ]
    in_format = "tsv"
    out_format = "tsv"
    CSVTK_SPLIT_TSV_TSV ( input, in_format, out_format )
}

workflow test_csvtk_split_csv_csv {

    input = [
        [ id:'test' ], // meta map
        [ file("https://raw.githubusercontent.com/SusiJo/test-datasets/modules/data/generic/csv/test.csv", checkIfExists: true) ]
    ]
    in_format = "csv"
    out_format = "csv"
    CSVTK_SPLIT_CSV_CSV ( input, in_format, out_format )
}
