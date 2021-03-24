#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KB_REF } from '../../../../software/kb/ref/main.nf' addParams( options: [:] )

workflow test_kb_ref {
    
    KB_REF ( 
        file("${launchDir}/tests/data/gencode.v26.chr21.GRCh38.p10.genome.fa", checkIfExists: true),
        file("${launchDir}/tests/data/gencode.v26.chr21.annotation.gtf", checkIfExists: true) 
        )
}
