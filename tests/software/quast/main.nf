#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QUAST }   from '../../../software/quast/main.nf'   addParams(options: [:])

workflow test_quast_ref {
    consensus = file("${launchDir}/tests/data/fasta/sarscov2/GCA_011545545.1_ASM1154554v1_cds_from_genomic.fna", checkIfExists: true)
    gff = file("${launchDir}/tests/data/gff/sarscov2/GCA_011545545.1_ASM1154554v1_genomic.gtf", checkIfExists: true)
    fasta = file("${launchDir}/tests/data/fasta/sarscov2/GCA_011545545.1_ASM1154554v1_genomic.fna", checkIfExists: true)
    def use_fasta = true
    def use_gtf = true

    QUAST( consensus, fasta, use_fasta, gff, use_gtf)
}

workflow test_quast_noref {
    consensus = file("${launchDir}/tests/data/fasta/sarscov2/GCA_011545545.1_ASM1154554v1_cds_from_genomic.fna", checkIfExists: true)
    gff = file('gff_dummy')
    fasta = file('fasta_dummy')
    def use_fasta = false
    def use_gtf = false

    QUAST( consensus, fasta, use_fasta, gff, use_gtf)
}
