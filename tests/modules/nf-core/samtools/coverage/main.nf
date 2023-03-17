#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_COVERAGE } from '../../../../../modules/nf-core/samtools/coverage/main.nf'
include { SAMTOOLS_SORT     } from '../../../../../modules/nf-core/samtools/sort/main.nf'

workflow test_samtools_coverage {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]
    SAMTOOLS_SORT ( input )
    SAMTOOLS_COVERAGE ( SAMTOOLS_SORT.out.bam )
}
