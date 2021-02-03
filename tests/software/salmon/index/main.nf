#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SALMON_INDEX } from '../../../software/salmon/index/main.nf' addParams( options: [args: '--minAssignedFrags 1'] )

workflow test_salmon_index {
    def genome_fasta     = file("${launchDir}/tests/data/fasta/sarscov2/GCA_011545545.1_ASM1154554v1_genomic.fna", checkIfExists: true)
    def transcript_fasta = file("${launchDir}/tests/data/fasta/sarscov2/GCA_011545545.1_ASM1154554v1_cds_from_genomic.fna", checkIfExists: true)
    SALMON_INDEX ( genome_fasta, transcript_fasta )
}