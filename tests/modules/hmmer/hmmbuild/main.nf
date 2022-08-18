#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMMER_HMMBUILD } from '../../../../modules/hmmer/hmmbuild/main.nf'

workflow test_hmmer_hmmbuild {
    
    input = [
        [ id: 'PF14720' ],      // meta map
        file('https://raw.githubusercontent.com/nf-core/test-datasets/phyloplace/testdata/PF14720_seed.alnfaa', checkIfExists: true)
    ]

    HMMER_HMMBUILD ( input )
}
