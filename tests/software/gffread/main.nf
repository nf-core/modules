#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GFFREAD } from '../../../software/gffread/main.nf' addParams( options: [:] )

workflow test_gffread {
    def input = []
    input = [ file("${launchDir}/tests/data/genomics/sarscov2/gtf/MN908947.3.gff3", checkIfExists: true) ]
    GFFREAD ( input )
}
