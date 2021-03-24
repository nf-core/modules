#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_ANNOTATEPEAKS } from '../../../../software/homer/annotatepeaks/main.nf' addParams( options: [:] )

workflow test_homer_annotatepeaks {

    def fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
    def gtf = file("${launchDir}/tests/data/genomics/sarscov2/gtf/test_genome.gtf", checkIfExists: true)

    def input = []
    input = [ [ id:'test'],
                file("${launchDir}/tests/data/genomics/sarscov2/bed/test.bed", checkIfExists: true) ]



    HOMER_ANNOTATEPEAKS ( input, fasta, gtf)
}
