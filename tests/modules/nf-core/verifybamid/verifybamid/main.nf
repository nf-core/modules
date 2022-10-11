#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VERIFYBAMID_VERIFYBAMID } from '../../../../../modules/nf-core/verifybamid/verifybamid/main.nf'

workflow test_verifybamid_verifybamid {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]

    refvcf = file(params.test_data['homo_sapiens']['illumina']['genmod_vcf_gz'], checkIfExists: true)

    VERIFYBAMID_VERIFYBAMID ( input, refvcf )
}
