#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FAIDX } from '../../../../software/samtools/faidx/main.nf' addParams( options: [:] )

workflow test_samtools_faidx {

    SAMTOOLS_FAIDX ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genomic.fasta", checkIfExists: true) )
}
