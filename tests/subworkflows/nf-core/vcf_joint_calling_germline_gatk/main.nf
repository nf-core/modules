#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_JOINT_CALLING_GERMLINE_GATK } from '../../../../subworkflows/nf-core/vcf_joint_calling_germline_gatk/main.nf'

workflow test_vcf_joint_calling_germline_gatk {

    // gvcf             // channel: [ val(meta), [ gvcf ], [ gvcf_index ], intervals]
    // fasta            // channel: [ val(meta), /path/to/reference/fasta]
    // fai              // channel: [ val(meta), /path/to/reference/fasta/index]
    // dict             // channel: [ val(meta), /path/to/reference/fasta/dict]
    // dbsnp            // channel: [ val(meta), /path/to/dbsnp/vcf]
    // dbsnp_tbi        // channel: [ val(meta), /path/to/dbsnp/vcf/index]

    input = [
        [ 'id': test],
        [ file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true)],
        [ file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'], checkIfExists: true)],
        file(params.test_data['homo_sapiens']['genome']['genome_interval_list'], checkIfExists: true)
    ]

    fasta = [ [ id:'genome', single_end:false ], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]
    fai   = [ [ id:'genome', single_end:false ], file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true) ]
    dict  = [ [ id:'genome', single_end:false ], file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true) ]

    dbsnp     = [ [ id:'dbsnp' ], file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true) ]
    dbsnp_tbi = [ [ id:'dbsnp' ], file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true) ]

    VCF_JOINT_CALLING_GERMLINE_GATK ( input, fasta, fai, dict, dbsnp, dbsnp_tbi )
}
