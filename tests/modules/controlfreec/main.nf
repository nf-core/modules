#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CONTROLFREEC_SOMATIC } from '../../../modules/controlfreec/main.nf'

workflow test_controlfreec {

    input = [
        [ id:'test', single_end:false, sex:'XX' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_mpileup'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_mpileup'], checkIfExists: true),
        [],[],[],[]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)

    dbsnp = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz'], checkIfExists: true)
    dbsnp_tbi = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz_tbi'], checkIfExists: true)

    CONTROLFREEC_SOMATIC (  input,
                            fasta,
                            fai,
                            [],
                            dbsnp,
                            dbsnp_tbi,
                            "/Users/monarchy/Projects/Coding/test-datasets/data/genomics/homo_sapiens/genome/chr21/sequence/genome.sizes",
                            "/Users/monarchy/Projects/Coding/test-datasets/data/genomics/homo_sapiens/genome/chr21/sequence/",
                            [],
                            [])
}
