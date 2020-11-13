#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { BEDTOOLS_GENOMECOV } from '../main.nf' addParams( options: [publish_dir:'test_bed_file'] )
// Define input channels

// Run the workflow
workflow test_bed_file {
    def input = []
    input = [ [ id:'test'], 
              [ file("${baseDir}/input/JK2067_downsampled_s0.1.bam", checkIfExists: true),] ]

    BEDTOOLS_GENOMECOV ( 
    input,
    file("${baseDir}/input/genome.sizes", checkIfExists: true) )

}

workflow {
    test_bed_file()
}


