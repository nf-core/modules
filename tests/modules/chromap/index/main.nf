#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHROMAP_INDEX } from '../../../../modules/chromap/index/main.nf' addParams( options: [:] )

workflow test_chromap_index {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    CHROMAP_INDEX ( input )
}
