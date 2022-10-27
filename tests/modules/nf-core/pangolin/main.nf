#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { PANGOLIN } from "$moduleDir/modules/nf-core/pangolin/main.nf"

workflow test_pangolin {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ] ]

    PANGOLIN ( input )
}
