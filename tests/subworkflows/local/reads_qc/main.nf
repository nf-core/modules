#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { READS_QC } from '../../../../subworkflows/local/reads_qc/main'

workflow test_reads_qc {

    trimmed_reads = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_umi_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_umi_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    untrimmed_reads = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'],     checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_2_fastq_gz'],     checkIfExists: true)
        ]
    ]

    READS_QC ( untrimmed_reads, trimmed_reads )

}
