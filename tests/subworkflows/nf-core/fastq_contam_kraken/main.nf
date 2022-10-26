#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQ_CONTAM_KRAKEN } from '../../../../subworkflows/nf-core/fastq_contam_kraken/main.nf'

workflow test_fastq_contam_kraken {

    input = [
        [ id:'test', single_end:true ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
    ]

    db = file("https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20220908.tar.gz")

    FASTQ_CONTAM_KRAKEN ( input , 25000, db)
}
