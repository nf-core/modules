#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LOFREQ_CALL } from '../../../../modules/lofreq/call/main.nf'

workflow test_lofreq_call {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    LOFREQ_CALL ( input, fasta )
}
