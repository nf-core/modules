#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LISSERO } from '../../../modules/lissero/main.nf'

workflow test_lissero {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file("https://github.com/MDU-PHL/LisSero/raw/master/tests/test_seq/NC_002973.fna", checkIfExists: true) ]

    LISSERO ( input )
}
