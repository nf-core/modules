#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC as FASTQC_SE } from '../main.nf'  addParams( options: [ publish_dir:'test_single_end' ] )
include { FASTQC as FASTQC_PE } from '../main.nf'  addParams( options: [ publish_dir:'test_paired_end' ] )

/*
 * Test with single-end data
 */
workflow test_single_end {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${baseDir}/input/test_single_end.fastq.gz", checkIfExists: true) ] ]

    FASTQC_SE ( input )
}

/*
 * Test with paired-end data
 */
workflow test_paired_end {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${baseDir}/input/test_R1.fastq.gz", checkIfExists: true),
                file("${baseDir}/input/test_R2.fastq.gz", checkIfExists: true) ] ]

    FASTQC_PE ( input )
}

workflow {
    test_single_end()
    test_paired_end()
}
