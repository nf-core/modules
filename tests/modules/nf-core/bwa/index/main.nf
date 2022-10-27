#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { BWA_INDEX } from "$moduleDir/modules/nf-core/bwa/index/main.nf"

workflow test_bwa_index {
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BWA_INDEX ( fasta )
}
