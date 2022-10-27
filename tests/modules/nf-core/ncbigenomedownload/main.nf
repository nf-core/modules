#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { NCBIGENOMEDOWNLOAD } from "$moduleDir/modules/nf-core/ncbigenomedownload/main.nf"

workflow test_ncbigenomedownload {

    input = [ [ id:'test', single_end:false ] ]

    accessions = []

    NCBIGENOMEDOWNLOAD ( input, accessions)
}


