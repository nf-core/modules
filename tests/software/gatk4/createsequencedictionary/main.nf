#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../../software/gatk4/createsequencedictionary/main.nf' addParams( options: [:] )

workflow test_gatk4_createsequencedictionary {
    GATK4_CREATESEQUENCEDICTIONARY ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) )
}
