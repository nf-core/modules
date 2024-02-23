#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE_BUILD } from '../../../../../modules/nf-core/bowtie/build/main.nf'

workflow test_bowtie_build {

    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BOWTIE_BUILD ( fasta )
}
