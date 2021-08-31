#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_BBMAPINDEX } from '../../../../modules/bbmap/bbmapindex/main.nf' addParams( options: [:] )

workflow test_bbmap_bbmapindex {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    BBMAP_BBMAPINDEX ( input )
}
