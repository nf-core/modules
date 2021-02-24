#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FAIDX } from '../../../../software/samtools/faidx/main.nf' addParams( options: [:] )

workflow test_samtools_faidx {

    SAMTOOLS_FAIDX ( file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true) )
}
