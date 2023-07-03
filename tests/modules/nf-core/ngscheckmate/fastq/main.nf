#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NGSCHECKMATE_FASTQ } from '../../../../../modules/nf-core/ngscheckmate/fastq/main.nf'

workflow test_ngscheckmate_fastq {

    input = [
        [ id:'test', single_end:true ], // meta map
            [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
    ]

    snp_pt = file("https://raw.githubusercontent.com/parklab/NGSCheckMate/master/SNP/SNP.pt", checkIfExists: true)

    NGSCHECKMATE_FASTQ ( input, snp_pt )
}

workflow test_ngscheckmate_fastq_paired {

    input = [
        [ id:'test', single_end:false ], // meta map
            [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ],
            [ file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
    ]

    snp_pt = file("https://raw.githubusercontent.com/parklab/NGSCheckMate/master/SNP/SNP.pt", checkIfExists: true)

    NGSCHECKMATE_FASTQ ( input, snp_pt )
}
