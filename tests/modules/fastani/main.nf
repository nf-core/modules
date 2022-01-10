#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTANI } from '../../../modules/fastani/main.nf'

workflow test_fastani {

    query = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    reference = file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)

    FASTANI ( [[ id:'test' ], query], reference )
}
