#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NCBIGENOMEDOWNLOAD } from '../../../../modules/nf-core/ncbigenomedownload/main.nf'

workflow test_ncbigenomedownload {

    input = [ [ id:'test', single_end:false ] ]

    accessions = []
    taxids = []
    groups = 'bacteria'

    NCBIGENOMEDOWNLOAD ( input, accessions, taxids, groups )
}


