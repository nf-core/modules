#!/usr/bin/env nextflow



include { BOWTIE_BUILD } from '../../../../modules/bowtie/build/main.nf'

workflow test_bowtie_build {
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BOWTIE_BUILD ( fasta )
}
