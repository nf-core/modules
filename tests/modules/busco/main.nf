#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO as BUSCO_BACTE } from '../../../modules/busco/main.nf'

// This tests genome decompression, empty input channels and data download
workflow test_busco {
   input = [ [ id:'test' ], file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true) ]    
   BUSCO_BACTE ( input,
                 "genome",
                 [],
                 [] )
}

