#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOMEPREPARATION } from '../../../../software/bismark/genomepreparation/main.nf' addParams( options: [:] )

workflow test_bismark_genomepreparation {

    BISMARK_GENOMEPREPARATION ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) )
}
