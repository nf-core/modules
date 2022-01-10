#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from '../../../../modules/bowtie2/build/main.nf'
include { BOWTIE2_ALIGN } from '../../../../modules/bowtie2/align/main.nf'

workflow test_bowtie2_align_single_end {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    save_unaligned = false

    BOWTIE2_BUILD ( fasta )
    BOWTIE2_ALIGN ( input, BOWTIE2_BUILD.out.index, save_unaligned )
}

workflow test_bowtie2_align_paired_end {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    save_unaligned = false

    BOWTIE2_BUILD ( fasta )
    BOWTIE2_ALIGN ( input, BOWTIE2_BUILD.out.index, save_unaligned )
}
