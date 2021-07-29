#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MALT_RUN } from '../../../../modules/malt/run/main.nf' addParams( options: [:] )

workflow test_malt_run {

    input = file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    fasta = file(params.test_data['sarscov2']['illumina']['genome_fasta'], checkIfExists: true)

    MALT_BUILD( fastas )
    MALT_RUN ( input, MALT_BUILD.out.index )
}

workflow test_malt_run {

    input = file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    fasta = file(params.test_data['sarscov2']['illumina']['genome_fasta'], checkIfExists: true)

    MALT_RUN ( input, MALT_BUILD.out.index )
}
