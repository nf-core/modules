#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_NORM } from '../../../../modules/bcftools/norm/main.nf' addParams( options: [:] )

workflow test_bcftools_norm {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    BCFTOOLS_NORM ( input )
}
