#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KHMER_NORMALIZEBYMEDIAN } from '../../../../modules/khmer/normalizebymedian/main.nf' addParams( options: [:] )

workflow test_khmer_normalizebymedian {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    KHMER_NORMALIZEBYMEDIAN ( input )
}
