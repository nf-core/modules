#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_COMPOSESTRTABLEFILE } from '../../../../../modules/nf-core/gatk4/composestrtablefile/main.nf'

workflow test_gatk4_composestrtablefile {
    
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_COMPOSESTRTABLEFILE ( fasta, fasta_fai, dict )
}
