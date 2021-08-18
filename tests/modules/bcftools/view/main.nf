#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_VIEW } from '../../../../modules/bcftools/view/main.nf' addParams( options: [:] )

workflow test_bcftools_view {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    BCFTOOLS_VIEW ( input )
}
