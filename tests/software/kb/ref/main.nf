#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KB_REF } from '../../../../software/kb/ref/main.nf' addParams( options: [:] )

workflow test_kb_ref {
    
    KB_REF ( 
        file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true),
        file("${launchDir}/tests/data/genomics/sarscov2/gtf/test_genome.gtf", checkIfExists: true) 
        )
}
