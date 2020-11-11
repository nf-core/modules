#!/usr/bin/env nextflow

nextflow.preview.dsl = 2
// Test -l, -r provided.
include { BEDTOOLS_SLOP } from '../main.nf' addParams( options: [publish_dir:'test_bed_file'] )
// Define input channels

// Run the workflow
workflow test_bed_file {
    def input = []
    input = [ [ id:'test', l:1, r:10 ], 
              [ file("${baseDir}/input/A.bed", checkIfExists: true),] ]

    BEDTOOLS_SLOP(
        input, 
        file("${baseDir}/input/genome.sizes", checkIfExists: true)
    )

}

workflow {
    test_bed_file()
}


