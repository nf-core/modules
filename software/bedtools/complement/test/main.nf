#!/usr/bin/env nextflow

nextflow.preview.dsl = 2
include { BEDTOOLS_COMPLEMENT } from '../main.nf' addParams( options: [publish_dir:'test_bed_file'] )
// Define input channels

// Run the workflow
workflow test_bed_file {
    def input = []
    input = [ [ id:'test'], 
              [ file("${baseDir}/input/A.bed", checkIfExists: true),] ]

    BEDTOOLS_COMPLEMENT(
        input, 
        file("${baseDir}/input/genome.sizes", checkIfExists: true)
    )

}

workflow {
    test_bed_file()
}


