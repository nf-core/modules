#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQ_CREATE_UMI_CONSENSUS_FGBIO } from '../../../../subworkflows/nf-core/fastq_create_umi_consensus_fgbio/main.nf'

workflow test_fastq_create_umi_consensus_fgbio_single_umi {

    reads = [
        [ id:'test_single', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_umi_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_umi_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta          =    file(params.test_data['homo_sapiens']['genome']['genome_fasta'],            checkIfExists: true)
    dict           =    file(params.test_data['homo_sapiens']['genome']['genome_dict'] ,            checkIfExists: true)

    FASTQ_CREATE_UMI_CONSENSUS_FGBIO ( reads, fasta, false, dict, "Adjacency", "bwa-mem", false )
}



workflow test_fastq_create_umi_consensus_fgbio_duplex_umi {

    reads = [
        [ id:'test_duplex', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_duplex_umi_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_duplex_umi_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta          =    file(params.test_data['homo_sapiens']['genome']['genome_fasta'],            checkIfExists: true)
    dict           =    file(params.test_data['homo_sapiens']['genome']['genome_dict'] ,            checkIfExists: true)

    FASTQ_CREATE_UMI_CONSENSUS_FGBIO ( reads, fasta, false, dict, "paired", "bwa-mem", true )
}
