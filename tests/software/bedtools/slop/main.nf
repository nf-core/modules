#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_SLOP } from '../../../../software/bedtools/slop/main.nf' addParams( options: [ suffix: '.slop', args: '-l 15 -r 30' ] )

workflow test_bedtools_slop {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true) ]

    BEDTOOLS_SLOP ( input, file("${launchDir}/tests/data/bed/genome.sizes", checkIfExists: true) )
}

