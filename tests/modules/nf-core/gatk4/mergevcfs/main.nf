#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MERGEVCFS } from '../../../../../modules/nf-core/gatk4/mergevcfs/main.nf'

workflow test_gatk4_mergevcfs {
    input = [ [ id:'test' ], // meta map
            [ file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz'], checkIfExists: true) ]
            ]

    dict  = [ [id:'dict'], file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)]

    GATK4_MERGEVCFS ( input, dict )
}

workflow test_gatk4_mergevcfs_no_dict {
    input = [ [ id:'test' ], // meta map
            [ file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz'], checkIfExists: true) ]
            ]

    GATK4_MERGEVCFS ( input, [[id:'genome'],[]] )
}

workflow test_gatk4_mergevcfs_no_dict_stubs {
    input = [ [ id:'test' ], // meta map
            [ file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz'], checkIfExists: true) ]
            ]

    GATK4_MERGEVCFS ( input, [[id:'genome'],[]] )
}
