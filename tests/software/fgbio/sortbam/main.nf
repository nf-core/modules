#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FGBIO_SORTBAM } from '../../../../software/fgbio/sortbam/main.nf' addParams( options: [:] )

workflow test_fgbio_sortbam {
    
    def input = []
    input = [ [ id:'test' ], // meta map
              file('https://github.com/fulcrumgenomics/fgbio/blob/master/src/test/resources/com/fulcrumgenomics/bam/200reads.bam?raw=true', checkIfExists: true) ]

    FGBIO_SORTBAM ( input )
}

