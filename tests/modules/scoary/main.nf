#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SCOARY } from '../../../modules/scoary/main.nf'

workflow test_scoary {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file("https://github.com/AdmiralenOla/Scoary/raw/master/scoary/exampledata/Gene_presence_absence.csv", checkIfExists: true),
              file("https://github.com/AdmiralenOla/Scoary/raw/master/scoary/exampledata/Tetracycline_resistance.csv", checkIfExists: true) ]

    tree = []
    SCOARY ( input, tree)
}
