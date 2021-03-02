#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_CREATESEQUENCEDICTIONARY } from '../../../../software/gatk/createsequencedictionary/main.nf' addParams( options: [:] )

workflow test_gatk_createsequencedictionary {
    GATK_CREATESEQUENCEDICTIONARY ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) )
}
