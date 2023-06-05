#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_FILTERVARIANTTRANCHES } from '../../../../../modules/nf-core/gatk4/filtervarianttranches/main.nf'
include { GATK4_CNNSCOREVARIANTS } from '../../../../../modules/nf-core/gatk4/cnnscorevariants/main.nf'
include { GATK4_HAPLOTYPECALLER  } from '../../../../../modules/nf-core/gatk4/haplotypecaller/main.nf'
workflow test_gatk4_filtervarianttranches {

    resources = [ file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true) ]
    resources_index  = [
                    file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true),
                ]

    input = [ [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz_tbi'], checkIfExists: true),
                []
            ]

    fasta  = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai    = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict   = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_FILTERVARIANTTRANCHES (input , resources, resources_index, fasta, fai, dict)
}
