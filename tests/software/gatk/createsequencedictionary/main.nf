#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_CREATESEQUENCEDICTIONARY } from '../../../../software/gatk/createsequencedictionary/main.nf' addParams( options: [:] )

workflow test_gatk_createsequencedictionary {
    GATK_CREATESEQUENCEDICTIONARY ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/GCA_011545545.1_ASM1154554v1_genomic.fasta", checkIfExists: true) )
}
