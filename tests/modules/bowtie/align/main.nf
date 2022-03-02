#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE_BUILD } from '../../../../modules/bowtie/build/main.nf'
include { BOWTIE_ALIGN } from '../../../../modules/bowtie/align/main.nf'

workflow test_bowtie_align_single_end {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BOWTIE_BUILD ( fasta )
    BOWTIE_ALIGN ( input, BOWTIE_BUILD.out.index )
}

workflow test_bowtie_align_paired_end {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BOWTIE_BUILD ( fasta )
    BOWTIE_ALIGN ( input, BOWTIE_BUILD.out.index )
}
