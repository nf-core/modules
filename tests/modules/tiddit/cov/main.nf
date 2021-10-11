#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TIDDIT_COV } from '../../../../modules/tiddit/cov/main.nf' addParams( options: [:] )

workflow test_tiddit_cov {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    TIDDIT_COV ( input )
}
