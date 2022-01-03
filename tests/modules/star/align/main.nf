#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE               } from '../../../../modules/star/genomegenerate/main.nf'
include { STAR_ALIGN                        } from '../../../../modules/star/align/main.nf'
include { STAR_ALIGN as STAR_FOR_ARRIBA     } from '../../../../modules/star/align/main.nf'
include { STAR_ALIGN as STAR_FOR_STARFUSION } from '../../../../modules/star/align/main.nf'

workflow test_star_alignment_single_end {
    input = [
        [ id:'test', single_end:true ], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_1_fastq_gz'], checkIfExists: true) ]
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf   = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    star_ignore_sjdbgtf = false
    seq_platform = 'illumina'
    seq_center = false

    STAR_GENOMEGENERATE ( fasta, gtf )
    STAR_ALIGN ( input, STAR_GENOMEGENERATE.out.index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
}

workflow test_star_alignment_paired_end {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf   = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    star_ignore_sjdbgtf = false
    seq_platform = 'illumina'
    seq_center = false

    STAR_GENOMEGENERATE ( fasta, gtf )
    STAR_ALIGN ( input, STAR_GENOMEGENERATE.out.index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
}


workflow test_star_alignment_paired_end_for_fusion {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf   = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    star_ignore_sjdbgtf = false
    seq_platform = 'illumina'
    seq_center = false

    STAR_GENOMEGENERATE ( fasta, gtf )
    STAR_FOR_ARRIBA ( input, STAR_GENOMEGENERATE.out.index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
}

workflow test_star_alignment_paired_end_for_starfusion {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf   = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    star_ignore_sjdbgtf = false
    seq_platform = false
    seq_center = false

    STAR_GENOMEGENERATE ( fasta, gtf )
    STAR_FOR_STARFUSION ( input, STAR_GENOMEGENERATE.out.index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
}
