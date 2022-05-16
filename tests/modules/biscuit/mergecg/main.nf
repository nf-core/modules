#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISCUIT_INDEX   } from '../../../../modules/biscuit/index/main.nf'
include { BISCUIT_MERGECG } from '../../../../modules/biscuit/mergecg/main.nf'

workflow test_biscuit_mergecg {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/biscuit/test-cg.bed.gz', checkIfExists: true)
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BISCUIT_INDEX( fasta )
    BISCUIT_MERGECG ( input, BISCUIT_INDEX.out.index )
}
