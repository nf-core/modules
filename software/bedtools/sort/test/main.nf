#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { BEDTOOLS_SORT } from '../main.nf' addParams( options: [publish_dir:'test_bed_file'] )
// Define input channels

// Run the workflow
workflow test_bed_file {
    def input = []
    input = [ [ id:'test'], 
              [ file("${baseDir}/input/A.bed", checkIfExists: true),] ]

    BEDTOOLS_SORT( input )

}

workflow {
    test_bed_file()
}


