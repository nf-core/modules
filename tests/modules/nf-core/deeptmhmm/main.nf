#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPTMHMM } from '../../../../modules/nf-core/deeptmhmm/main.nf'

workflow test_deeptmhmm {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    DEEPTMHMM ( fasta )
}

workflow test_deeptmhmm_gz {

    fasta_gz = [ file(params.test_data['sarscov2']['genome']['genome_fasta_gz'], checkIfExists: true) ]

    DEEPTMHMM ( fasta_gz )
}
