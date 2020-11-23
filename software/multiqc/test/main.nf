#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MULTIQC } from '../main.nf'  addParams( options: [ publish_dir:'test_single_end' ] )

/*
 * Test with single-end data
 */
workflow test_single_end {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${baseDir}/input/test_single_end.fastq.gz", checkIfExists: true) ] ]

    FASTQC_SE ( input )
}

workflow {
    test_single_end()
}
