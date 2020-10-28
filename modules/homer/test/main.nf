#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_CONFIGUREHOMER } from '../configurehomer/main.nf'  addParams( options: [ publish_dir:'test_single_end' ] )
include { HOMER_MAKETAGDIRECTORY } from '../maketagdirectory/main.nf' addParams( options: [ publish_dir:'test_one_file' ] )
include { HOMER_MAKETAGDIRECTORY as HOMER_MAKETAGDIRECTORY_TWO } from '../maketagdirectory/main.nf' addParams( options: [ publish_dir:'test_two_file' ] )
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_TWO } from '../annotatepeaks/main.nf' addParams( options: [ publish_dir:'test_two_file' ] )

/*
 * Test with single-end data
 */
workflow test_one_file {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${baseDir}/input/A.bed", checkIfExists: true) ] ]

    HOMER_CONFIGUREHOMER()
    HOMER_MAKETAGDIRECTORY( input )
}

workflow test_two_file {

    def input2 = []
    input2 = [ [ id:'test_two', single_end:true ], // meta map
              [ file("${baseDir}/input/A.bed", checkIfExists: true),
                file("${baseDir}/input/B.bed", checkIfExists: true) ] ]

    // HOMER_CONFIGUREHOMER( )
    HOMER_MAKETAGDIRECTORY_TWO( input2 )
    HOMER_FINDPEAKS( HOMER_MAKETAGDIRECTORY_TWO.out.tagdir )
    // FIXME HOMER_ANNOTATEPEAKS_TWO( HOMER_MAKETAGDIRECTORY_TWO.out.tagdir )

}

workflow {
    // FIXME test_one_file()
    test_two_file()
}
