#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_CONFIGUREHOMER } from '../configurehomer/main.nf'  addParams( options: [ publish_dir:'test_single_end' ] )
include { HOMER_MAKETAGDIRECTORY } from '../maketagdirectory/main.nf' addParams( options: [ publish_dir:'test_single_end' ] )

/*
 * Test with single-end data
 */
workflow test {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${baseDir}/input/A.bed", checkIfExists: true) ] ]

    HOMER_CONFIGUREHOMER()
    HOMER_MAKETAGDIRECTORY( input )
}

workflow {
    test()
}
