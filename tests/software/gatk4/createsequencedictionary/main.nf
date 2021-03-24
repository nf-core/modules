#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../../software/gatk4/createsequencedictionary/main.nf' addParams( options: [:] )

workflow test_gatk4_createsequencedictionary {
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)

    GATK4_CREATESEQUENCEDICTIONARY ( fasta )
}
