#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QUAST }   from '../../../software/quast/main.nf'   addParams(options: [:])

workflow test_quast_ref {
    consensus = file("${launchDir}/tests/data/genomics/sarscov2/fasta/GCA_011545545.1_ASM1154554v1_cds_from_genomic.fna", checkIfExists: true)
    gff = file("${launchDir}/tests/data/genomics/sarscov2/gtf/GCA_011545545.1_ASM1154554v1_genomic.gtf", checkIfExists: true)
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/GCA_011545545.1_ASM1154554v1_genomic.fasta", checkIfExists: true)
    def use_fasta = true
    def use_gtf = true

    QUAST( consensus, fasta, gff, use_fasta, use_gtf )
}

workflow test_quast_noref {
    consensus = file("${launchDir}/tests/data/genomics/sarscov2/fasta/GCA_011545545.1_ASM1154554v1_cds_from_genomic.fna", checkIfExists: true)
    gff = file('gff_dummy')
    fasta = file('fasta_dummy')
    def use_fasta = false
    def use_gtf = false

    QUAST( consensus, fasta, gff, use_fasta, use_gtf )
}
