#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP } from '../../../software/gunzip/main.nf' addParams( options: [:] )

workflow test_gunzip {

    def input = []
    input = [ file("${launchDir}/tests/data/vcf/test.vcf.gz", checkIfExists: true) ]

    GUNZIP ( input )
}
