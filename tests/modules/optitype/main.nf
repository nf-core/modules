#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { OPTITYPE } from '../../../modules/optitype/main.nf'

workflow test_optitype {
    input = [ [ id:'test', seq_type:'dna' ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_hla_sorted_bam'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_hla_sorted_bam_bai'], checkIfExists: true)
            ]

    OPTITYPE ( input )
}
