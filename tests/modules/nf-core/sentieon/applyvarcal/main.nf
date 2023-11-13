#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SENTIEON_APPLYVARCAL } from '../../../../../modules/nf-core/sentieon/applyvarcal/main.nf'

workflow test_sentieon_applyvarcal {
    input = [ [ id:'test'], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test2_recal'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test2_recal_idx'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test2_tranches'], checkIfExists: true)
            ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)

    SENTIEON_APPLYVARCAL( input, [ [:], fasta], [ [:], fai] )
}
