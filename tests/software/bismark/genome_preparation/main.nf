#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOME_PREPARATION } from '../../../../software/bismark/genome_preparation/main.nf' addParams( options: [:] )

workflow test_bismark_genome_preparation {

    BISMARK_GENOME_PREPARATION ( file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true) )
}
