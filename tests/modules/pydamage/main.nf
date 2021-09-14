#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYDAMAGE } from '../../../modules/pydamage/main.nf' addParams( options: [:] )

workflow test_pydamage {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]

    PYDAMAGE ( input )
}
