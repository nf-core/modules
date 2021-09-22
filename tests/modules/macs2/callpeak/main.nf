#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MACS2_CALLPEAK } from '../../../../modules/macs2/callpeak/main.nf' addParams( options: ['args': '--bdg', ] )
include { MACS2_CALLPEAK as MACS2_CALLPEAK_BROAD } from '../../../../modules/macs2/callpeak/main.nf' addParams( options: ['args': '--bdg --broad', ] )

workflow test_macs2_callpeak {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]
    gsize = 29903

    MACS2_CALLPEAK ( input, gsize )
}

workflow test_macs2_callpeak_broad {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]
    gsize = 29903

    MACS2_CALLPEAK_BROAD ( input, gsize )
}
