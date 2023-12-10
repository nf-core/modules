#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC } from '../../../fastqc/main'
include { MULTIQC } from "../../main"

workflow test_multiqc {
    input = [
        [ id: 'test', single_end: false ],
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)]
    ]

    FASTQC  ( input )
    MULTIQC ( FASTQC.out.zip.collect { it[1] }, [], [], [] )
}
