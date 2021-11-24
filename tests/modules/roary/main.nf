#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ROARY } from '../../../modules/roary/main.nf' addParams( options: [:] )

workflow test_roary {
    
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("https://github.com/bactopia/bactopia-tests/raw/main/data/species/portiera/genome/gff/GCF_000292685.gff", checkIfExists: true),
                file("https://github.com/bactopia/bactopia-tests/raw/main/data/species/portiera/genome/gff/GCF_000298385.gff", checkIfExists: true),
                file("https://github.com/bactopia/bactopia-tests/raw/main/data/species/portiera/genome/gff/GCF_002849995.gff", checkIfExists: true) ]
    ]

    ROARY ( input )
}
