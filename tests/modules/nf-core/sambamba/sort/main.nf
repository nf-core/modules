#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMBAMBA_SORT } from '../../../../../modules/nf-core/sambamba/sort/main.nf'

// test for bam output
workflow test_sambamba_sort{
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SAMBAMBA_SORT ( input )
}
