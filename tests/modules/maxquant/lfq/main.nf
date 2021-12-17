#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAXQUANT_LFQ } from '../../../../modules/maxquant/lfq/main.nf' addParams( options: [:] )

workflow test_maxquant_lfq {
        
    input = [ [ id:'test' ], // meta map
              file(params.fasta, checkIfExists: true), file(params.paramfile, checkIfExists: true) 
            ]


    rawfiles = [file('https://github.com/wombat-p/MaxQuant-Workflow/raw/dev/Nextflow/data_test/OVEMB150205_12.raw') , file('https://github.com/wombat-p/MaxQuant-Workflow/raw/dev/Nextflow/data_test/OVEMB150205_14.raw') ]

    MAXQUANT_LFQ ( input, rawfiles.collect() )
}
