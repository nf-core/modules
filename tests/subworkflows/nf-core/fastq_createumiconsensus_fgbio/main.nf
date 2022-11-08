#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQ_CREATEUMICONSENSUS_FGBIO } from '../../../../subworkflows/nf-core/fastq_createumiconsensus_fgbio/main.nf'

workflow test_fastq_createumiconsensus_fgbio_single_umi {

    reads = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_umi_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_umi_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta          =    file(params.test_data['homo_sapiens']['genome']['genome_fasta'],            checkIfExists: true)
    dict           =    file(params.test_data['homo_sapiens']['genome']['genome_dict'] ,            checkIfExists: true)
    read_structure =    "+T 12M11S+T"

    FASTQ_CREATEUMICONSENSUS_FGBIO ( reads, fasta, dict, read_structure, "Adjacency", "bwa-mem", false )
}



workflow test_fastq_createumiconsensus_fgbio_duplex_umi {

    reads = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_duplex_umi_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_duplex_umi_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta          =    file(params.test_data['homo_sapiens']['genome']['genome_fasta'],            checkIfExists: true)
    dict           =    file(params.test_data['homo_sapiens']['genome']['genome_dict'] ,            checkIfExists: true)
    read_structure =    "38S2M1S+T 2M34S+T"

    FASTQ_CREATEUMICONSENSUS_FGBIO ( reads, fasta, dict, read_structure, "paired", "bwa-mem", true )
}
