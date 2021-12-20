#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMBLASTER } from '../../../modules/samblaster/main.nf'

workflow test_samblaster {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_umi_unsorted_bam'], checkIfExists: true) ]

    SAMBLASTER ( input )
}
