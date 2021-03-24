#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UCSC_BED12TOBIGBED } from '../../../../software/ucsc/bed12tobigbed/main.nf' addParams( options: [:] )

workflow test_ucsc_bed12tobigbed {
    input = [ [ id: 'test' ], // meta map 
              [ file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed12", checkIfExists: true ) ] 
            ]
    sizes = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.sizes", checkIfExists: true)

    UCSC_BED12TOBIGBED ( input, sizes )
}
