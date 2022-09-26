#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRESEQ_LCEXTRAP } from '../../../../modules/preseq/lcextrap/main.nf'

//
// Test with single-end data
//
workflow test_preseq_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file('https://github.com/smithlabcode/preseq/raw/master/data/SRR1003759_5M_subset.mr', checkIfExists: true) ]
            ]
    PRESEQ_LCEXTRAP ( input )
}

//
// Test with paired-end data
//
workflow test_preseq_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file('https://github.com/smithlabcode/preseq/raw/master/data/SRR1003759_5M_subset.mr', checkIfExists: true) ]
            ]

    PRESEQ_LCEXTRAP ( input )
}

