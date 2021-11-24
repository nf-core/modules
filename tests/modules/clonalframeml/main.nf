#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CLONALFRAMEML } from '../../../modules/clonalframeml/main.nf' addParams( options: [:] )

workflow test_clonalframeml {
    
    input = [ [ id:'test' ], // meta map
              file("https://github.com/bactopia/bactopia-tests/raw/main/data/species/haemophilus_influenzae/genome/genome.aln.nwk", checkIfExists: true),
              file("https://github.com/bactopia/bactopia-tests/raw/main/data/species/haemophilus_influenzae/genome/genome.aln.gz", checkIfExists: true),]

    CLONALFRAMEML ( input )
}
