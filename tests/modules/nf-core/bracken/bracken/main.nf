#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR           } from '../../../../../modules/nf-core/untar/main.nf'
include { KRAKEN2_KRAKEN2 } from '../../../../../modules/nf-core/kraken2/kraken2/main.nf'
include { BRACKEN_BRACKEN } from '../../../../../modules/nf-core/bracken/bracken/main.nf'

workflow test_bracken_bracken_single_end_default_args {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    db    = file(params.test_data['sarscov2']['genome']['kraken2_bracken_tar_gz'], checkIfExists: true)

    ch_db = UNTAR ( [[:], db] ).untar
        .map { it[1] }
    KRAKEN2_KRAKEN2 ( input, ch_db, false, false )
    BRACKEN_BRACKEN ( KRAKEN2_KRAKEN2.out.report, ch_db )
}

workflow test_bracken_bracken_paired_end_default_args {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    db    = file(params.test_data['sarscov2']['genome']['kraken2_bracken_tar_gz'], checkIfExists: true)

    ch_db = UNTAR ( [[:], db] ).untar
        .map { it[1] }
    KRAKEN2_KRAKEN2 ( input, ch_db, false, false )
    BRACKEN_BRACKEN ( KRAKEN2_KRAKEN2.out.report, ch_db )
}
