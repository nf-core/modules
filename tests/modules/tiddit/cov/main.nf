#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TIDDIT_COV } from '../../../../modules/tiddit/cov/main.nf'

workflow test_tiddit_cov {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    TIDDIT_COV ( input, fasta )
}

workflow test_tiddit_cov_no_ref {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]

    TIDDIT_COV ( input, [] )
}
