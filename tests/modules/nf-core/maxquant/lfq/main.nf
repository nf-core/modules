#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAXQUANT_LFQ } from '../../../../../modules/nf-core/maxquant/lfq/main.nf'

workflow test_maxquant_lfq {
        
    input = [ [ id:'test' ], // meta map
              file(params.test_data['proteomics']['database']['yeast_ups'], checkIfExists: true), file(params.test_data['proteomics']['parameter']['maxquant'] , checkIfExists: true)
            ]


    rawfiles = [file(params.test_data['proteomics']['msspectra']['ups_file1']) , file(params.test_data['proteomics']['msspectra']['ups_file2'])]

    MAXQUANT_LFQ ( input, rawfiles.collect() )
}
