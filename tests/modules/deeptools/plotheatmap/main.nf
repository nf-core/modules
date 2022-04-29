#!/usr/bin/env nextflow



include { DEEPTOOLS_PLOTHEATMAP } from '../../../../modules/deeptools/plotheatmap/main.nf'

workflow test_deeptools_plotheatmap {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_computematrix_mat_gz'], checkIfExists: true)
            ]

    DEEPTOOLS_PLOTHEATMAP ( input )
}
