#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_BBMAPINDEX } from '../../../../modules/bbmap/bbmapindex/main.nf' addParams( options: [:] )

workflow test_bbmap_bbmapindex {
    
    input = [ [ id:'test' ], // meta map
              file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]

    BBMAP_BBMAPINDEX ( input )
}
