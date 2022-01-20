#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MASH_DIST } from '../../../../modules/mash/dist/main.nf'

workflow test_mash_dist {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]

    reference =  file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    MASH_DIST ( input, reference )
}
