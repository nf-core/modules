#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HISAT2_BUILD } from '../../../../software/hisat2/build/main.nf' addParams( options: [:] )

workflow test_hisat2_build {
    
    def input = []
    input = [ file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true),
            file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.gtf", checkIfExists: true),
            file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
        ]

    HISAT2_BUILD ( input )
}
