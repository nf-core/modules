#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NCBIGENOMEDOWNLOAD } from '../../../modules/ncbigenomedownload/main.nf' addParams( options: [ args: '-A GCF_000013425.1 --formats genbank,fasta,assembly-stats bacteria '] )

workflow test_ncbigenomedownload {
    
    input = [ [ id:'test', single_end:false ] ]

    accessions = []

    NCBIGENOMEDOWNLOAD ( input, accessions)
}


