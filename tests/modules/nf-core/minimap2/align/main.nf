#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINIMAP2_ALIGN } from '../../../../../modules/nf-core/minimap2/align/main.nf'

workflow test_minimap2_align_single_end {
    input = [
        [ id:'test', single_end:true ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    bam_format = true
    cigar_paf_format = false
    cigar_bam = false

    MINIMAP2_ALIGN ( input, fasta, bam_format, cigar_paf_format, cigar_bam)
}

workflow test_minimap2_align_paired_end {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    bam_format  = true
    cigar_paf_format = false
    cigar_bam = false

    MINIMAP2_ALIGN ( input, fasta, bam_format, cigar_paf_format, cigar_bam )
}

workflow test_minimap2_align_self_align {
    input = [
        [ id:'test', single_end:true ], // meta map
        [ file(params.test_data['sarscov2']['genome']['genome_fasta_gz'], checkIfExists: true) ]
    ]
    fasta = []
    bam_format = true
    cigar_paf_format = false
    cigar_bam = false

    MINIMAP2_ALIGN ( input, fasta, bam_format, cigar_paf_format, cigar_bam)
}
