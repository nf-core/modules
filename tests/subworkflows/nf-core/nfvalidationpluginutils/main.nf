#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NFVALIDATIONPLUGINUTILS } from '../../../../subworkflows/nf-core//main.nf'

workflow test_nfvalidationpluginutils {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    NFVALIDATIONPLUGINUTILS ( input )
}
