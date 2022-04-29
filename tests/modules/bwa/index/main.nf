#!/usr/bin/env nextflow



include { BWA_INDEX } from '../../../../modules/bwa/index/main.nf'

workflow test_bwa_index {
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BWA_INDEX ( fasta )
}
