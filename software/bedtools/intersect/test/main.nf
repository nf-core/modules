#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { BEDTOOLS_INTERSECT } from '../main.nf' addParams( options: [publish_dir:'test_bed_file'] )
// Run the workflow
workflow test_bed_file {
    def input = []
    input = [ [ id:'test'], 
              [ file("${baseDir}/input/A.bed", checkIfExists: true),
                file("${baseDir}/input/B.bed", checkIfExists: true) ] ]

    BEDTOOLS_INTERSECT( input )

}

workflow {
    test_bed_file()
}


