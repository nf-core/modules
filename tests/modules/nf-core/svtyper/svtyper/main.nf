#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SVTYPER_SVTYPER } from '../../../../../modules/nf-core/svtyper/svtyper/main.nf'
include { GRIDSS_GRIDSS   } from '../../../../../modules/nf-core/gridss/gridss/main.nf'
include { BWA_INDEX       } from '../../../../../modules/nf-core/bwa/index/main.nf'

workflow test_svtyper_svtyper {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        []
    ]
    fasta = [[id:'fasta'],file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[id:'fasta_fai'],file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]

    GRIDSS_GRIDSS( input, fasta, fasta_fai, BWA_INDEX(fasta).index )

    bam = [
        [ id:'bam' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
    ]

    SVTYPER_SVTYPER ( bam, GRIDSS_GRIDSS.out.vcf, fasta, fasta_fai )
}
