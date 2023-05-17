#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRINITY } from '../../../../modules/nf-core/trinity/main.nf'

workflow test_trinity_paired_end {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_1_fastq_gz'], checkIfExists: true),
          file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_2_fastq_gz'], checkIfExists: true) ]
    ]

    TRINITY ( input )
}


workflow test_trinity_single_end {
    
    input = [
        [ id:'test', single_end:true ], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_1_fastq_gz'], checkIfExists: true) ]
    ]

    TRINITY ( input )
}