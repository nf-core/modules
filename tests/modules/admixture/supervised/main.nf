#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ADMIXTURE_SUPERVISED } from '../../../../modules/admixture/supervised/main.nf' addParams( options: [:] )

workflow test_admixture_supervised {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    ADMIXTURE_SUPERVISED ( input )
}
