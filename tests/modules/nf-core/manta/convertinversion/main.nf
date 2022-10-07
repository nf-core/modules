#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MANTA_CONVERTINVERSION } from '../../../../../modules/nf-core/manta/convertinversion/main.nf'
include { MANTA_TUMORONLY        } from '../../../../../modules/nf-core/manta/tumoronly/main.nf'

workflow test_manta_convertinversion {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram_crai'], checkIfExists: true),
        [], []
    ]

    fasta   = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai     = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)

    MANTA_TUMORONLY ( input, fasta, fai )

    MANTA_CONVERTINVERSION ( MANTA_TUMORONLY.out.tumor_sv_vcf, fasta )
}
