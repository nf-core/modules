#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_SENDSKETCH } from '../../../../../modules/nf-core/bbmap/sendsketch/main.nf'

workflow test_bbmap_sendsketch {
    
    input = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    BBMAP_SENDSKETCH ( input )
}
