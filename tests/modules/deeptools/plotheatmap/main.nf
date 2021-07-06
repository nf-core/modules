#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPTOOLS_PLOTHEATMAP } from '../../../../software/deeptools/plotheatmap/main.nf' addParams( options: [:] )

workflow test_deeptools_plotheatmap {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_computematrix_mat_gz'], checkIfExists: true)
            ]

    DEEPTOOLS_PLOTHEATMAP ( input )
}
