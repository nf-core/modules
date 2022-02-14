#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_COMBINEGVCFS } from '../../../../modules/gatk4/combinegvcfs/main.nf'

workflow test_gatk4_combinegvcfs {
    
    input = [
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test2_vcf'], checkIfExists: true)
    ]
    
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    GATK4_COMBINEGVCFS ( fasta, input )
}

