#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_GERMLINECNVCALLER } from '../../../../modules/gatk4/germlinecnvcaller/main.nf'

workflow test_gatk4_germlinecnvcaller {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    // model = file('/fs1/resources/ref/hg38/gatk_cnv/cnvref/ploidy-model/', checkIfExists: true)

    GATK4_COLLECTREADCOUNTS ( input, fasta, fai, dict )
    GATK4_DETERMINEGERMLINECONTIGPLOIDY ( GATK4_COLLECTREADCOUNTS.out.tsv)//, model )
    GATK4_GERMLINECNVCALLER ( GATK4_COLLECTREADCOUNTS.out.tsv, GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.tar )
}
