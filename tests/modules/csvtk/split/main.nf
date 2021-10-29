#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CSVTK_SPLIT as CSVTK_SPLIT_TSV_TSV } from '../../../../modules/csvtk/split/main.nf' addParams( options: [args: "-C '&' --fields '#chr' "])
include { CSVTK_SPLIT as CSVTK_SPLIT_CSV_CSV } from '../../../../modules/csvtk/split/main.nf' addParams( options: [args: "-C '&' --fields '#chr' "])

workflow test_csvtk_split_tsv_tsv {

    input = [
        [ id:'test' ], // meta map
        [ file("https://raw.githubusercontent.com/nf-core/test-datasets/epitopeprediction/testdata/variants/variants.tsv", checkIfExists: true) ]
    ]
    in_format = "tsv"
    out_format = "tsv"
    CSVTK_SPLIT_TSV_TSV ( input, in_format, out_format )
}

workflow test_csvtk_split_csv_csv {

    input = [
        [ id:'test' ], // meta map
        [ file("https://raw.githubusercontent.com/SusiJo/test-datasets/modules/data/generic/txt/variants.csv", checkIfExists: true) ]
    ]
    in_format = "csv"
    out_format = "csv"
    CSVTK_SPLIT_CSV_CSV ( input, in_format, out_format )
}
