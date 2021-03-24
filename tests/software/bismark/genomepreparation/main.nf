#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOMEPREPARATION } from '../../../../software/bismark/genomepreparation/main.nf' addParams( options: [:] )

workflow test_bismark_genomepreparation {
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)

    BISMARK_GENOMEPREPARATION ( fasta )
}
