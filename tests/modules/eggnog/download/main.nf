#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EGGNOG_DOWNLOAD } from '../../../../modules/eggnog/download/main.nf' addParams( options: [:] )

workflow test_eggnog_download {
    
    input = [ [ id:'test'] ]

    EGGNOG_DOWNLOAD ( input )
}
