#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_ANNOTATEPEAKS } from '../../../../software/homer/annotatepeaks/main.nf' addParams( options: [:] )

workflow test_homer_annotatepeaks {
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed", checkIfExists: true) 
            ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    gtf   = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.gtf", checkIfExists: true)

    HOMER_ANNOTATEPEAKS ( input, fasta, gtf )
}
