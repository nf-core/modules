#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FAIDX } from '../../../../software/samtools/faidx/main.nf' addParams( options: [:] )

workflow test_samtools_faidx {

    SAMTOOLS_FAIDX ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/GCA_011545545.1_ASM1154554v1_genomic.fasta", checkIfExists: true) )
}
