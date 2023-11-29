#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GAMMA_GAMMA } from '../../../../../modules/nf-core/gamma/gamma/main.nf'

workflow test_unzip {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true),
    ]

    db = [ file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/srst2/ResGANNCBI_20210507_srst2.fasta", checkIfExists: true), ]

    GAMMA_GAMMA ( input, db )
}

workflow test_gamma {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    db = [ file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true) ]

    GAMMA_GAMMA ( input, db )
}
