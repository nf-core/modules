#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRESTO_FILTERSEQ } from '../../../../../modules/nf-core/presto/filterseq/main.nf'

workflow test_presto_filterseq {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['fastq']['test_airrseq_1_fastq_gz']}, checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['fastq']['test_airrseq_2_fastq_gz']}, checkIfExists: true)
    ]

    PRESTO_FILTERSEQ ( input )
}
