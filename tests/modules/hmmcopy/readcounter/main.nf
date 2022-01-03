#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMMCOPY_READCOUNTER } from '../../../../modules/hmmcopy/readcounter/main.nf'

workflow test_hmmcopy_readcounter {

    input = [ [ id:'test'], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)]
            ]
    HMMCOPY_READCOUNTER ( input )
}
