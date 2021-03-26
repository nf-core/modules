#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../../software/gatk4/createsequencedictionary/main.nf' addParams( options: [:] )

workflow test_gatk4_createsequencedictionary {
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    GATK4_CREATESEQUENCEDICTIONARY ( fasta )
}
