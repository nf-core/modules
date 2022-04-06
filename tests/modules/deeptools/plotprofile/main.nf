#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPTOOLS_PLOTPROFILE } from '../../../../modules/deeptools/plotprofile/main.nf'

workflow test_deeptools_plotprofile {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_computematrix_mat_gz'], checkIfExists: true)
            ]

    DEEPTOOLS_PLOTPROFILE ( input )
}
