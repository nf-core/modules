#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { PRODIGAL } from "$moduleDir/modules/nf-core/prodigal/main.nf"

workflow test_prodigal {
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
            ]

    PRODIGAL ( input , "gff")
}
