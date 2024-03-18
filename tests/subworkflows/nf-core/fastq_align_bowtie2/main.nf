#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from '../../../../modules/nf-core/bowtie2/build/main.nf'
include { FASTQ_ALIGN_BOWTIE2 } from '../../../../subworkflows/nf-core/fastq_align_bowtie2/main.nf'


workflow test_align_bowtie2_single_end {
    input = [
        [ id:'test', single_end:true ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    save_unaligned = false
    sort = false

    BOWTIE2_BUILD ( fasta )
    FASTQ_ALIGN_BOWTIE2 ( input, BOWTIE2_BUILD.out.index, save_unaligned, sort, [[],[]] )
}

workflow test_align_bowtie2_paired_end {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    save_unaligned = false
    sort = false

    BOWTIE2_BUILD ( fasta )
    FASTQ_ALIGN_BOWTIE2 ( input, BOWTIE2_BUILD.out.index, save_unaligned, sort, [[],[]] )
}

