#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DRAGEN } from '../../../../modules/nf-core/dragen/main.nf'

workflow test_dragen {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ],
    ]
    reference = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    DRAGEN ( input, reference )
}
