#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_CREATESOMATICPANELOFNORMALS } from '../../../../modules/gatk4/createsomaticpanelofnormals/main.nf' addParams( options: [:] )

workflow test_gatk4_createsomaticpanelofnormals {

    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true) , file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz'], checkIfExists: true)] ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    GATK4_CREATESOMATICPANELOFNORMALS ( input , fasta)
}
