#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_SLOP } from '../../../../software/bedtools/slop/main.nf' addParams( options: [args: '-l 15 -r 30'] )

workflow test_bedtools_slop {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/bed/sarscov2.bed", checkIfExists: true) ]

    BEDTOOLS_SLOP ( input, file("${launchDir}/tests/data/genomics/sarscov2/bed/sarscov2_genome.sizes", checkIfExists: true) )
}

