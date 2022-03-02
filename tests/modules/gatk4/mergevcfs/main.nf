#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MERGEVCFS } from '../../../../modules/gatk4/mergevcfs/main.nf'

workflow test_gatk4_mergevcfs {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test2_vcf'], checkIfExists: true) ] 
            ]
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_MERGEVCFS ( input, dict, false )
}

workflow test_gatk4_mergevcfs_refdict {
    def input = []
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test2_vcf'], checkIfExists: true) ] 
            ]
    dict  = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_MERGEVCFS ( input, dict, true )
}
