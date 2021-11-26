#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_INDEX } from '../../../../modules/bbmap/index/main.nf'

workflow test_bbmap_index {
    
    input = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    BBMAP_INDEX ( input )
}
