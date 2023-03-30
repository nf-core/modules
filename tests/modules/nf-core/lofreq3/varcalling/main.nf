#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LOFREQ3_VARCALLING } from '../../../../../modules/nf-core/lofreq3/varcalling/main.nf'

workflow test_lofreq3_varcalling {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    LOFREQ3_VARCALLING ( input, false )
}
