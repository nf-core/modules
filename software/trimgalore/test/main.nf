#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRIMGALORE } from '../main.nf'

/*
 * Test with single-end data
 */
workflow test_single_end {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${baseDir}/input/test_single_end.fastq.gz", checkIfExists: true) ] ]

    TRIMGALORE ( input, [ publish_dir:'test_single_end' ] )
}

/*
 * Test with paired-end data
 */
workflow test_paired_end {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${baseDir}/input/test_R1.fastq.gz", checkIfExists: true),
                file("${baseDir}/input/test_R2.fastq.gz", checkIfExists: true) ] ]

    TRIMGALORE ( input, [ publish_dir:'test_paired_end' ] )
}

workflow {
    test_single_end()
    test_paired_end()
}
