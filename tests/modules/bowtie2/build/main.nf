#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from '../../../../modules/bowtie2/build/main.nf'

workflow test_bowtie2_build {
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BOWTIE2_BUILD ( fasta )
}
