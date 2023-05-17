#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR           } from '../../../../../modules/nf-core/untar/main.nf'
include { KRAKEN2_KRAKEN2 } from '../../../../../modules/nf-core/kraken2/kraken2/main.nf'
include { BRACKEN_BRACKEN } from '../../../../../modules/nf-core/bracken/bracken/main.nf'
include { BRACKEN_COMBINEBRACKENOUTPUTS } from '../../../../../modules/nf-core/bracken/combinebrackenoutputs/main.nf'

workflow test_bracken_combinebrackenoutputs {

    input = Channel.of(
        [[ id:'test', single_end:false, threshold:0, taxonomic_level:'G', read_length:100 ], [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true), file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]],
        [[ id:'test2', single_end:false, threshold:0, taxonomic_level:'G', read_length:100 ], [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true), file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]],
    )
    db    = file(params.test_data['sarscov2']['genome']['kraken2_bracken_tar_gz'], checkIfExists: true)

    ch_db = UNTAR ( [[:], db] ).untar
        .map { it[1] }
    KRAKEN2_KRAKEN2 ( input, ch_db, false, false )
    BRACKEN_BRACKEN ( KRAKEN2_KRAKEN2.out.report, ch_db )

    ch_input_for_combinebrackenouputs = BRACKEN_BRACKEN.out.reports
                                            .map{it[1]}
                                            .collect()
                                            .map{ [ [id: 'db'], it ] }

    BRACKEN_COMBINEBRACKENOUTPUTS ( ch_input_for_combinebrackenouputs )
}

