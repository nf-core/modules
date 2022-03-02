#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METHYLDACKEL_MBIAS } from '../../../../modules/methyldackel/mbias/main.nf'

workflow test_methyldackel_mbias {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam_bai'], checkIfExists: true) ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    METHYLDACKEL_MBIAS ( input, fasta, fai )
}
