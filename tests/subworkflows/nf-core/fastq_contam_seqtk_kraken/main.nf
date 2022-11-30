#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQ_CONTAM_SEQTK_KRAKEN } from '../../../../subworkflows/nf-core/fastq_contam_seqtk_kraken/main.nf'
include { UNTAR                     } from '../../../modules/nf-core/untar/main.nf'

workflow test_FASTQ_CONTAM_SEQTK_KRAKEN {

    input = [
        [ id:'test', single_end:true ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ],
    ]

    db    = [ [], file(params.test_data['sarscov2']['genome']['kraken2_tar_gz'], checkIfExists: true) ]

    UNTAR ( db )

    FASTQ_CONTAM_SEQTK_KRAKEN ( Channel.of(input) , Channel.of(25000), UNTAR.out.untar.map{ it[1] })
}

workflow test_FASTQ_CONTAM_SEQTK_KRAKEN_MULTIPLE_N {

    input = [
        [ id:'test', single_end:true ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ],
    ]

    db    = [ [], file(params.test_data['sarscov2']['genome']['kraken2_tar_gz'], checkIfExists: true) ]

    UNTAR ( db )

    FASTQ_CONTAM_SEQTK_KRAKEN ( Channel.of(input) , Channel.of(12500, 25000), UNTAR.out.untar.map{ it[1] })
}
