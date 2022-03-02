#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MLST } from '../../../modules/mlst/main.nf'

workflow test_mlst {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file("https://raw.githubusercontent.com/nf-core/test-datasets/bactmap/genome/NCTC13799.fna", checkIfExists: true) ]

    MLST ( input )
}
