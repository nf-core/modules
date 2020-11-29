#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MULTIQC } from '../main.nf'  addParams( options: [ publish_dir:'test_multi' ] )

/*
 * Test with single-end data
 */
workflow test_multi {

    def input = []
    input = [
        file("${baseDir}/input/test_1_fastqc.html", checkIfExists: true),
        file("${baseDir}/input/test_2_fastqc.html", checkIfExists: true),
        file("${baseDir}/input/test_1_fastqc.zip", checkIfExists: true),
        file("${baseDir}/input/test_2_fastqc.zip", checkIfExists: true)
    ]

    MULTIQC ( input )
}

workflow {
    test_multi()
}
