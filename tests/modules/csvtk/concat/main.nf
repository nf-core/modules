#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CSVTK_CONCAT } from '../../../../modules/csvtk/concat/main.nf'

workflow test_csvtk_concat {

    input = [ 
        [ id:'test' ], // meta map
        [ file("https://github.com/nf-core/test-datasets/raw/bacass/bacass_hybrid.csv", checkIfExists: true),
          file("https://github.com/nf-core/test-datasets/raw/bacass/bacass_long.csv", checkIfExists: true),
          file("https://github.com/nf-core/test-datasets/raw/bacass/bacass_short.csv", checkIfExists: true) ]
    ]

    in_format = "tsv"
    out_format = "csv"

    CSVTK_CONCAT ( input, in_format, out_format )
}
