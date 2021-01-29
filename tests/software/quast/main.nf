#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QUAST as QUAST_NOREF } from '../../../software/quast/quast_noref/main.nf' addParams(options: [:])
include { QUAST as QUAST_REF }   from '../../../software/quast/quast_ref/main.nf'   addParams(options: [:])

workflow test_quast_noref {
    consensus = file("${launchDir}/tests/data/fasta/sarscov2/GCA_011545545.1_ASM1154554v1_genomic.fna", checkIfExists: true)

    QUAST_NOREF( consensus )
}

workflow test_quast_ref {
    consensus = file("${launchDir}/tests/data/fasta/sarscov2/GCA_011545545.1_ASM1154554v1_cds_from_genomic.fna", checkIfExists: true)
    gff = file("${launchDir}/tests/data/gff/sarscov2/GCA_011545545.1_ASM1154554v1_genomic.gtf", checkIfExists: true)
    fasta = file("${launchDir}/tests/data/fasta/sarscov2/GCA_011545545.1_ASM1154554v1_genomic.fna", checkIfExists: true)

    QUAST_REF( consensus, fasta, gff )
}
