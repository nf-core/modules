#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

include BEDTOOLS_SLOPREFSEQ from '../sloprefseq/main.nf' addParams( options: [publish_dir:'test_bed_file'])

// Define input channels

// Run the workflow
workflow test_bed_file{
    def input = []
    input = [ [ id:'test' ], 
              [ file("${baseDir}/input/A.bed", checkIfExists: true),
                file("${baseDir}/input/genome.sizes", checkIfExists: true)] ]

    BEDTOOLS_SLOPREFSEQ ( input )

}

workflow {
    test_bed_file()
}
