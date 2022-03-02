#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LOFREQ_CALLPARALLEL } from '../../../../modules/lofreq/callparallel/main.nf'

workflow test_lofreq_callparallel {
    
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    LOFREQ_CALLPARALLEL ( input, fasta, fai )
}
