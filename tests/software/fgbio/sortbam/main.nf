#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FGBIO_SORTBAM } from '../../../../software/fgbio/sortbam/main.nf' addParams( options: [:] )

workflow test_fgbio_sortbam {
    
    def input = []
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['bam']['test_paired_end.sorted.bam'], checkIfExists: true) ]
    FGBIO_SORTBAM ( input )
}

