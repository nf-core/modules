#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HICAP } from '../../../modules/hicap/main.nf' addParams( options: [:] )

workflow test_hicap {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file("https://github.com/bactopia/bactopia-tests/raw/main/data/species-specific/haemophilus-influenzae/GCF_900478275.fna.gz", checkIfExists: true) ]
    
    database_dir = []
    model_fp = []

    HICAP ( input, database_dir, model_fp )
}
