#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE_BUILD } from '../../../../software/bowtie/build/main.nf' addParams( options: [:] )

workflow test_bowtie_build {
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
    BOWTIE_BUILD ( fasta )
}
