#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARIANTBAM } from '../../../modules/variantbam/main.nf'

workflow test_variantbam {

    input = [ [ id:'test', single_end:false ],
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]

    VARIANTBAM ( input )
}
