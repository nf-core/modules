#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE } from '../../../../modules/nf-core/star/genomegenerate/main.nf'
include { FASTQ_ALIGN_STAR    } from '../../../../subworkflows/nf-core/fastq_align_star/main.nf'


workflow test_align_star_single_end {
    input = [
        [ id:'test', single_end:true ], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
    ]
    fasta = [
        [ id:'test_fasta', single_end:true ], // meta map
        [ file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]
    ]
    gtf = [
        [ id:'test_gtf', single_end:true ], // meta map
        [ file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true) ]
    ]
    star_ignore_sjdbgtf = true
    seq_platform = 'illumina'
    seq_center = false

    STAR_GENOMEGENERATE ( fasta, gtf )
    FASTQ_ALIGN_STAR ( input, STAR_GENOMEGENERATE.out.index, [], star_ignore_sjdbgtf, seq_platform, seq_center, [[],[]] )
}

workflow test_align_star_paired_end {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [ id:'test_fasta', single_end:true ], // meta map
        [ file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]
    ]
    gtf = [
        [ id:'test_gtf', single_end:true ], // meta map
        [ file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true) ]
    ]
    star_ignore_sjdbgtf = true
    seq_platform = 'illumina'
    seq_center = false

    STAR_GENOMEGENERATE ( fasta, gtf )
    FASTQ_ALIGN_STAR ( input, STAR_GENOMEGENERATE.out.index, [], star_ignore_sjdbgtf, seq_platform, seq_center, [[],[]] )
}

