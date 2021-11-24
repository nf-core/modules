#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PIRATE } from '../../../modules/pirate/main.nf' addParams( options: [:] )

workflow test_pirate {
    
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("https://github.com/bactopia/bactopia-tests/raw/main/data/species/portiera/genome/gff/GCF_000292685.gff", checkIfExists: true),
                file("https://github.com/bactopia/bactopia-tests/raw/main/data/species/portiera/genome/gff/GCF_000298385.gff", checkIfExists: true),
                file("https://github.com/bactopia/bactopia-tests/raw/main/data/species/portiera/genome/gff/GCF_002849995.gff", checkIfExists: true) ]
    ]

    PIRATE ( input )
}
