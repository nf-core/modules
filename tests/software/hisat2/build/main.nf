#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HISAT2_BUILD } from '../../../../software/hisat2/build/main.nf' addParams( options: [:] )

workflow test_hisat2_build {
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    gtf = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.gtf", checkIfExists: true)
    splice = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.splice_sites.txt", checkIfExists: true)

    HISAT2_BUILD ( fasta, gtf, splice )
}
