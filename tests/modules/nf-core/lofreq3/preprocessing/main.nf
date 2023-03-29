#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LOFREQ3_PREPROCESSING } from '../../../../../modules/nf-core/lofreq3/preprocessing/main.nf'

workflow test_lofreq3 {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
    ]

    LOFREQ3_PREPROCESSING ( input )
}
