#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_SORT_SAMTOOLS as BAM_SORT_SAMTOOLS_SINGLE_END } from '../../../../subworkflows/nf-core/bam_sort_samtools/main'
include { BAM_SORT_SAMTOOLS as BAM_SORT_SAMTOOLS_PAIRED_END } from '../../../../subworkflows/nf-core/bam_sort_samtools/main'

workflow test_bam_sort_samtools_single_end {
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BAM_SORT_SAMTOOLS_SINGLE_END ( input, fasta )
}

workflow test_bam_sort_samtools_paired_end {
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BAM_SORT_SAMTOOLS_PAIRED_END ( input, fasta )
}
