#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQUENCES_NEWICK_PHYLOPLACE_EPANG_GAPPA } from '../../../../subworkflows/nf-core/sequences_newick_phyloplace_epang_gappa/main.nf'

workflow test_sequences_newick_phyloplace_epang_gappa {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    SEQUENCES_NEWICK_PHYLOPLACE_EPANG_GAPPA ( input )
}
