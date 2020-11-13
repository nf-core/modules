#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { BEDTOOLS_MERGE } from '../main.nf' addParams( options: [publish_dir:'test_bed_file'] )
// Define input channels
// TODO ensure input file is presorted (could use output of sort module)
// Run the workflow
workflow test_bed_file {
    def input = []
    input = [ [ id:'test'], 
              [ file("${baseDir}/input/A.bed", checkIfExists: true),] ]

    BEDTOOLS_MERGE( input )

}

workflow {
    test_bed_file()
}


