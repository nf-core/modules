#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR } from '../../../../../modules/nf-core/untar/main'
include { KRAKEN2_KRAKEN2 } from '../../../../../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_COMBINEKREPORTS } from '../../../../../modules/nf-core/krakentools/combinekreports/main.nf'

workflow test_krakentools_combinekreports {

    input = Channel.of(
        [[ id:'test', single_end:false ], [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true), file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]],
        [[ id:'test2', single_end:false ], [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true), file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]],
    )


    db    = [ [], file(params.test_data['sarscov2']['genome']['kraken2_tar_gz'], checkIfExists: true) ]

    UNTAR ( db )
    KRAKEN2_KRAKEN2 ( input, UNTAR.out.untar.map{ it[1] }, false, false )

    KRAKENTOOLS_COMBINEKREPORTS ( KRAKEN2_KRAKEN2.out.report.map{ [[id:"test"], it[1]] }.groupTuple())
}
