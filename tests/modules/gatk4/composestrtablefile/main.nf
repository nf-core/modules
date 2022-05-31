#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_COMPOSESTRTABLEFILE } from '../../../../modules/gatk4/composestrtablefile/main.nf'

workflow test_gatk4_composestrtablefile {
    
    input = [
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    ]

    GATK4_COMPOSESTRTABLEFILE ( input )
}
