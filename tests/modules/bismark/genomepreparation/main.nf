#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOMEPREPARATION } from '../../../../modules/bismark/genomepreparation/main.nf'

workflow test_bismark_genomepreparation {
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BISMARK_GENOMEPREPARATION ( fasta )
}
