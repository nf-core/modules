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

    chrfiles = file("/Users/monarchy/Projects/Coding/modules/sequences/", checkIfExists: true)
    target_bed = file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)

    CONTROLFREEC_SOMATIC (  input,
                            fasta,
                            fai,
                            [],
                            dbsnp,
                            dbsnp_tbi,
                            chrfiles,
                            [],
                            target_bed,
                            []
                        )
}
