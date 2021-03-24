#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_COMPLEMENT } from '../../../../software/bedtools/complement/main.nf' addParams( options: [suffix: '_out'] )

workflow test_bedtools_complement {
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed", checkIfExists: true)
            ]
    sizes = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.sizes", checkIfExists: true)

    BEDTOOLS_COMPLEMENT ( input, sizes )
}
