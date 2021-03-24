#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_SLOP } from '../../../../software/bedtools/slop/main.nf' addParams( options: [args: '-l 15 -r 30', suffix: '_out'] )

workflow test_bedtools_slop {
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed", checkIfExists: true)
            ]
    sizes = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.sizes", checkIfExists: true)

    BEDTOOLS_SLOP ( input, sizes )
}

