#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX } from '../../../../software/bwa/index/main.nf' addParams( options: [:] )

workflow test_bwa_index {
    BWA_INDEX ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genomic.fasta", checkIfExists: true) )
}
