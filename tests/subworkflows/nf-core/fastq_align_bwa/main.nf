#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX       } from '../../../../modules/nf-core/bwa/index/main.nf'
include { FASTQ_ALIGN_BWA } from '../../../../subworkflows/nf-core/fastq_align_bwa/main.nf'


workflow test_fastq_align_bwa_single_end {
    input = [
        [ id:'test', single_end:true ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    sort = false

    BWA_INDEX ( fasta )
    FASTQ_ALIGN_BWA ( input, BWA_INDEX.out.index, sort, [ ] )
}

workflow test_fastq_align_bwa_paired_end {
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
    sort = false

    BWA_INDEX ( fasta )
    FASTQ_ALIGN_BWA ( input, BWA_INDEX.out.index, sort, [ ] )
}
