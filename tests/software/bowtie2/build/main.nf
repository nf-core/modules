#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from '../../../software/bowtie2/build/main.nf' addParams( options: [:] )

workflow test_bowtie2_build {
    fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)
    BOWTIE2_BUILD ( fasta )
}