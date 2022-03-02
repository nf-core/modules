#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FGBIO_SORTBAM } from '../../../../modules/fgbio/sortbam/main.nf'

workflow test_fgbio_sortbam {
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
    FGBIO_SORTBAM ( input )
}
