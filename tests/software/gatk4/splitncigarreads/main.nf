#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_SPLITNCIGARREADS } from '../../../../software/gatk4/splitncigarreads/main.nf' addParams( options: [:] )

workflow test_gatk4_splitncigarreads {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/bam/sarscov2_aln.bam", checkIfExists: true)] ]

    fasta = file("tests/data/fasta/sarscov2/GCA_011545545.1_ASM1154554v1_genomic.fna", checkIfExists: true)
    fai = file("tests/data/fasta/sarscov2/GCA_011545545.1_ASM1154554v1_genomic.fna.fai", checkIfExists: true)
    dict = file("tests/data/fasta/sarscov2/GCA_011545545.1_ASM1154554v1_genomic.dict", checkIfExists: true)

    GATK4_SPLITNCIGARREADS ( input, fasta, fai, dict )
}
