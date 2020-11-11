#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

include { BEDTOOLS_SLOP } from '../slop/main.nf' addParams( options: [publish_dir:'test_bed_file'])
include { BEDTOOLS_REMOVEGENES } from '../removegenes/main.nf' addParams( options: [publish_dir:'test_bed_file'])

// Define input channels

// Run the workflow
workflow test_bed_file {
    def input = []
    input = [ [ id:'test', single_end:true ], 
              [ file("${baseDir}/input/A.bed", checkIfExists: true),] ]

    BEDTOOLS_SLOP(
        input, 
        file("${baseDir}/input/genome.sizes", checkIfExists: true)
    )
    BEDTOOLS_REMOVEGENES(
        BEDTOOLS_SLOP.out.slopbed,
        file("${baseDir}/input/B.metatranscripts", checkIfExists: true)
        )

}

workflow {
    test_bed_file()
}
