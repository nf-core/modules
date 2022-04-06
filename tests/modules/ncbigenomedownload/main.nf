#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NCBIGENOMEDOWNLOAD } from '../../../modules/ncbigenomedownload/main.nf'

workflow test_ncbigenomedownload {
    
    input = [ [ id:'test', single_end:false ] ]

    accessions = []

    NCBIGENOMEDOWNLOAD ( input, accessions)
}


