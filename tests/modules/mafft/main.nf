#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAFFT } from '../../../modules/mafft/main.nf'

workflow test_mafft {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['scaffolds_fasta'], checkIfExists: true)
    ]

    MAFFT ( input )
}
