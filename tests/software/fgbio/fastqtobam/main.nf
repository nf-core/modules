#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FGBIO_FASTQTOBAM } from '../../../../software/fgbio/fastqtobam/main.nf' addParams( options: [:] )

workflow test_fgbio_fastqtobam {
    
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.bam", checkIfExists: true) ]

    FGBIO_FASTQTOBAM ( input )
}
