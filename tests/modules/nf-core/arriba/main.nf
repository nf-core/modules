#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE   } from '../../../../modules/nf-core/star/genomegenerate/main.nf'
include { STAR_ALIGN            } from '../../../../modules/nf-core/star/align/main.nf'
include { ARRIBA                } from '../../../../modules/nf-core/arriba/main.nf'

workflow test_arriba_single_end {

    input = [ [ id:'test', single_end:true ], // meta map
                [   file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_1_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = Channel.fromPath(params.test_data['homo_sapiens']['genome']['genome_fasta']).map { it -> [[id:it.Name], it] }.collect()
    fai = Channel.fromPath(params.test_data['homo_sapiens']['genome']['genome_fasta_fai']).map { it -> [[id:it.Name], it] }.collect()
    gtf = Channel.fromPath(params.test_data['homo_sapiens']['genome']['genome_gtf']).map { it -> [[id:it.Name], it] }.collect()
    star_ignore_sjdbgtf = false
    seq_platform = 'illumina'
    seq_center = false

    STAR_GENOMEGENERATE ( fasta, gtf )
    STAR_ALIGN ( input, STAR_GENOMEGENERATE.out.index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
    ARRIBA ( STAR_ALIGN.out.bam, fasta, gtf , [], [], [], [], [])
}

workflow test_arriba_paired_end {

    input = [ [ id:'test', single_end:false ], // meta map
                [   file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_1_fastq_gz'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_2_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = Channel.fromPath(params.test_data['homo_sapiens']['genome']['genome_fasta']).map { it -> [[id:it.Name], it] }.collect()
    fai = Channel.fromPath(params.test_data['homo_sapiens']['genome']['genome_fasta_fai']).map { it -> [[id:it.Name], it] }.collect()
    gtf = Channel.fromPath(params.test_data['homo_sapiens']['genome']['genome_gtf']).map { it -> [[id:it.Name], it] }.collect()
    star_ignore_sjdbgtf = false
    seq_platform = 'illumina'
    seq_center = false

    STAR_GENOMEGENERATE ( fasta, gtf )
    STAR_ALIGN ( input, STAR_GENOMEGENERATE.out.index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
    ARRIBA ( STAR_ALIGN.out.bam, fasta, gtf, [], [], [], [], [])
}
