#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAFFT              } from '../../../../modules/nf-core/mafft/main.nf'
include { MAFFT as MAFFT_ADD } from '../../../../modules/nf-core/mafft/main.nf'

workflow test_mafft {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['scaffolds_fasta'], checkIfExists: true)
    ]

    MAFFT ( input, [] )
}

workflow test_mafft_add {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/epang/reference.alnfaa.gz', checkIfExists: true),
    ]

    MAFFT_ADD (
        input,
        file('https://raw.githubusercontent.com/nf-core/test-datasets/phyloplace/testdata/PF14720_3_sequences.faa', checkIfExists: true)
    )
}
