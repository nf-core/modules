#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MERGEVCFS } from "$moduleDir/modules/nf-core/gatk4/mergevcfs/main.nf"

workflow test_gatk4_mergevcfs {
    def input = []
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test2_vcf'], checkIfExists: true) ]
            ]

    dict  = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_MERGEVCFS ( input, dict )
}

workflow test_gatk4_mergevcfs_no_dict {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test2_vcf'], checkIfExists: true) ]
            ]

    GATK4_MERGEVCFS ( input, [] )
}

workflow test_gatk4_mergevcfs_no_dict_stubs {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test2_vcf'], checkIfExists: true) ]
            ]

    GATK4_MERGEVCFS ( input, [] )
}
