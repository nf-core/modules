#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_CRAMSIZE } from '../../../../../modules/nf-core/samtools/cramsize/main.nf'

workflow test_samtools_cramsize {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true)
    ]

    SAMTOOLS_CRAMSIZE ( input )
}