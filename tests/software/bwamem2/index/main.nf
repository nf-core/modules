#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMEM2_INDEX } from '../../../../software/bwamem2/index/main.nf' addParams( options: [:] )

workflow test_bwamem2_index {
    BWAMEM2_INDEX ( file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true) )
}