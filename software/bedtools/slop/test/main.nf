#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BEDTOOLS_SLOP as BEDTOOLS_SLOP_S } from '../main.nf' addParams( options: [ publish_dir:'test_asymmetrical' ], header: true  )
include { BEDTOOLS_SLOP as BEDTOOLS_SLOP_AS } from '../main.nf' addParams( options: [ publish_dir:'test_symmetrical' ], header: true )

// To run with header and pct enabled, type --pct true and --header true with nextflow run command.

// Run the workflow
/*
Test with l/r method
*/
workflow test_asymmetrical {
    def input = []
    input = [ [ id:'test', symmetry: false ], 
              [ file("${baseDir}/input/A.bed", checkIfExists: true),] ]

    BEDTOOLS_SLOP_AS(
        input, 
        file("${baseDir}/input/genome.sizes", checkIfExists: true)
    )

}
/*
Test with b method
*/
workflow test_symmetrical {
    def input = []
    input = [ [ id:'test', symmetry: true], 
              [ file("${baseDir}/input/A.bed", checkIfExists: true),] ]

    BEDTOOLS_SLOP_S(
        input, 
        file("${baseDir}/input/genome.sizes", checkIfExists: true)
    )

}

workflow {
    test_asymmetrical()
    test_symmetrical()
}

// output has - b 10 and -l 1 -r 10 as test values. pct and header are both false.
