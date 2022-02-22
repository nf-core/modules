#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_COMBINEGVCFS } from '../../../../modules/gatk4/combinegvcfs/main.nf'

workflow test_gatk4_combinegvcfs {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
           ]

    
    vcffiles = [ file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true),
                 file(params.test_data['homo_sapiens']['illumina']['test2_genome_vcf'], checkIfExists: true) ]

    vcf_idx = [ file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_idx'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_genome_vcf_idx'], checkIfExists: true) ]  


//    vcffiles = [ file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true),
//                 file(params.test_data['homo_sapiens']['illumina']['test2_genome_vcf_gz'], checkIfExists: true) ]

    GATK4_COMBINEGVCFS ( input, vcffiles, vcf_idx )
}

