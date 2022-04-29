#!/usr/bin/env nextflow



include { SAMTOOLS_DEPTH } from '../../../../modules/samtools/depth/main.nf'

workflow test_samtools_depth {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true) ]

    SAMTOOLS_DEPTH ( input )
}
