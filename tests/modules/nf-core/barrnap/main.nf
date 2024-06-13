#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BARRNAP } from '../../../../modules/nf-core/barrnap/main.nf'

workflow test_barrnap {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]

    BARRNAP ( input )
}
