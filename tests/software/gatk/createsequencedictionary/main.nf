#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_CREATESEQUENCEDICTIONARY } from '../../../../software/gatk/createsequencedictionary/main.nf' addParams( options: [:] )

workflow test_gatk_createsequencedictionary {
    GATK_CREATESEQUENCEDICTIONARY ( file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true) )
}