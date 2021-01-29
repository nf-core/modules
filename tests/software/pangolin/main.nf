#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PANGOLIN } from '../../../software/pangolin/main.nf' addParams(options: [:])

workflow test_pangolin {
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/../GCA_011545545.1_ASM1154554v1_genomic.fna", checkIfExists: true) ] ]

    PANGOLIN( input )
}
