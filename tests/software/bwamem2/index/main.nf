#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMEM2_INDEX } from '../../../../software/bwamem2/index/main.nf' addParams( options: [:] )

workflow test_bwamem2_index {
    BWAMEM2_INDEX ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/GCA_011545545.1_ASM1154554v1_genomic.fasta", checkIfExists: true) )
}
