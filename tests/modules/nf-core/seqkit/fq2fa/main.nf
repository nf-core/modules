#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQKIT_FQ2FA } from '../../../../../modules/nf-core/seqkit/fq2fa/main.nf'

workflow test_seqkit_fq2fa {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
    ]

    SEQKIT_FQ2FA ( input )
}
