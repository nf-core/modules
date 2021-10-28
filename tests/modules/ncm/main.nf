#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NCM } from '../../../modules/ncm/main.nf' addParams( options: [:] )

workflow test_ncm {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    NCM ( input )
}
