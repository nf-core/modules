#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from '../../../../software/bowtie2/build/main.nf' addParams( options: [:] )

workflow test_bowtie2_build {
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/GCA_011545545.1_ASM1154554v1_genomic.fasta", checkIfExists: true)
    BOWTIE2_BUILD ( fasta )
}
