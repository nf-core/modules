#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_COMPLEMENT } from '../../../../software/bedtools/complement/main.nf' addParams( options: [suffix: '_out'] )

workflow test_bedtools_complement {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/bed/test.bed", checkIfExists: true) ]

    BEDTOOLS_COMPLEMENT ( input, file("${launchDir}/tests/data/genomics/sarscov2/bed/test.genome.sizes", checkIfExists: true) )
}
