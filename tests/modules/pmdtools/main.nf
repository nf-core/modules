#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PMDTOOLS } from '../../../modules/pmdtools/main.nf' addParams( options: [:] )

workflow test_pmdtools {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    PMDTOOLS ( input )
}
