#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX } from '../../../../software/bwa/index/main.nf' addParams( options: [:] )

workflow test_bwa_index {
    BWA_INDEX ( file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true) )
}