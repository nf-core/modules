#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_PREPROCESSINTERVALS } from '../../../../modules/gatk4/preprocessintervals/main.nf'

workflow test_gatk4_preprocessintervals {
    
    bed = file(params.test_data['homo_sapiens']['genome']['genome_blacklist_interval_bed'], checkIfExists: true)
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_PREPROCESSINTERVALS ( bed, fasta, fai, dict )
}
