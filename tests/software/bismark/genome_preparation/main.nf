#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOME_PREPARATION } from '../../../../software/bismark/genome_preparation/main.nf' addParams( options: [:] )

workflow test_bismark_genome_preparation {

    BISMARK_GENOME_PREPARATION ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genomic.fasta", checkIfExists: true) )
}
