#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MERGEMUTECTSTATS } from '../../../../modules/gatk4/mergemutectstats/main.nf'

workflow test_gatk4_mergemutectstats {

    input = [
        [ id:'test', single_end:false ], // meta map
         file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_vcf_gz_stats'], checkIfExists: true)
        ]

    GATK4_MERGEMUTECTSTATS ( input )
}
