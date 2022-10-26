#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQ_CONTAM_KRAKEN } from '../../../../subworkflows/nf-core/fastq_contam_kraken/main.nf'
include { UNTAR           }     from '../../../modules/nf-core/untar/main.nf'

workflow test_fastq_contam_kraken {

    input = [
        [ id:'test', single_end:true ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
    ]

    db    = [ [], file(params.test_data['sarscov2']['genome']['kraken2_tar_gz'], checkIfExists: true) ]

    UNTAR ( db )

    FASTQ_CONTAM_KRAKEN ( input , 25000, UNTAR.out.untar.map{ it[1] })
}
