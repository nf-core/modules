#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_COMPLEMENT } from '../../../../software/bedtools/complement/main.nf' addParams( options: [:] )


workflow test_bedtools_complement {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true),
              file("${launchDir}/tests/data/bed/genome.sizes", checkIfExists: true) ]

    BEDTOOLS_COMPLEMENT( input )
}

