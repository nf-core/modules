#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GFFREAD } from '../../../software/gffread/main.nf' addParams( options: [:] )

workflow test_gffread {
    input = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.gff3", checkIfExists: true)

    GFFREAD ( input )
}
