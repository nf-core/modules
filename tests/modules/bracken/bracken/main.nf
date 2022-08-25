#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR           } from '../../../../modules/untar/main.nf'
include { KRAKEN2_KRAKEN2 } from '../../../../modules/kraken2/kraken2/main.nf'
include { BRACKEN_BRACKEN } from '../../../../modules/bracken/bracken/main.nf'

workflow test_bracken_bracken_single_end_default_args {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    db    = file(params.test_data['sarscov2']['genome']['kraken2_bracken_tar_gz'], checkIfExists: true)

    ch_db = UNTAR ( [[:], db] ).untar
        .map { it[1] }
    KRAKEN2_KRAKEN2 ( input, ch_db )
    BRACKEN_BRACKEN ( KRAKEN2_KRAKEN2.out.txt, ch_db )
}

workflow test_bracken_bracken_single_end_custom_args {
    input = [ [ id:'test', single_end:true, threshold:0, taxonomic_level:'G', read_length:100 ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    db    = file(params.test_data['sarscov2']['genome']['kraken2_bracken_tar_gz'], checkIfExists: true)

    ch_db = UNTAR ( [[:], db] ).untar
        .map { it[1] }
    KRAKEN2_KRAKEN2 ( input, ch_db )
    BRACKEN_BRACKEN ( KRAKEN2_KRAKEN2.out.txt, ch_db )
}

workflow test_bracken_bracken_paired_end_default_args {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    db    = file(params.test_data['sarscov2']['genome']['kraken2_bracken_tar_gz'], checkIfExists: true)

    ch_db = UNTAR ( [[:], db] ).untar
        .map { it[1] }
    KRAKEN2_KRAKEN2 ( input, ch_db )
    BRACKEN_BRACKEN ( KRAKEN2_KRAKEN2.out.txt, ch_db )
}

workflow test_bracken_bracken_paired_end_custom_args {
    input = [ [ id:'test', single_end:false, threshold:0, taxonomic_level:'G', read_length:100 ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    db    = file(params.test_data['sarscov2']['genome']['kraken2_bracken_tar_gz'], checkIfExists: true)

    ch_db = UNTAR ( [[:], db] ).untar
        .map { it[1] }
    KRAKEN2_KRAKEN2 ( input, ch_db )
    BRACKEN_BRACKEN ( KRAKEN2_KRAKEN2.out.txt, ch_db )
}
