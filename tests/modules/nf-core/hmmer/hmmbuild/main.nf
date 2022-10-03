#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMMER_HMMBUILD } from '../../../../../modules/nf-core/hmmer/hmmbuild/main.nf'

workflow test_hmmer_hmmbuild {
    
    input = [
        [ id: 'PF14720' ],      // meta map
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/hmmer/PF14720_seed.alnfaa.gz', checkIfExists: true)
    ]

    HMMER_HMMBUILD ( input, [] )
}
