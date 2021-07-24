#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { JUPYTERNOTEBOOK } from '../../../modules/jupyternotebook/main.nf' addParams( options: [:] )

workflow test_jupyternotebook {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    JUPYTERNOTEBOOK ( input )
}
