#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SVABA         } from '../../../../modules/nf-core/svaba/main.nf'
include { BWA_INDEX     } from '../../../../modules/nf-core/bwa/index/main.nf'

workflow test_svaba_nocontrol {

    input     = [ [ id:'test' ], // meta map
                    file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                    [],
                    []
                ]
    fasta = [[id:'fasta'],file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[id:'fasta'],file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]

    BWA_INDEX (fasta)
    dbsnp     = [[id:'dbsnp'], file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true)]
    dbsnp_tbi = [[id:'dbsnp'],file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true)]

    SVABA ( input, fasta, fasta_fai, BWA_INDEX.out.index, dbsnp, dbsnp_tbi, [[],[]] )
}
workflow test_svaba {

    input     = [ [ id:'test' ], // meta map
                    file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)
                ]
    fasta = [[id:'fasta'],file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[id:'fasta'],file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]

    BWA_INDEX (fasta)
    dbsnp     = [[id:'dbsnp'], file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true)]
    dbsnp_tbi = [[id:'dbsnp'],file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true)]

    SVABA ( input, fasta, fasta_fai, BWA_INDEX.out.index, dbsnp, dbsnp_tbi, [[],[]] )
}
