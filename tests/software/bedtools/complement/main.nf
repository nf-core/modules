#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_COMPLEMENT } from '../../../../software/bedtools/complement/main.nf' addParams( options: [suffix:'.complement'] )

workflow test_bedtools_complement {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true) ]

    BEDTOOLS_COMPLEMENT ( input, file("${launchDir}/tests/data/bed/genome.sizes", checkIfExists: true) )
}
