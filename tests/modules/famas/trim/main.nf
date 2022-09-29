#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FAMAS_TRIM } from '../../../../modules/famas/trim/main.nf' addParams( options: [args: ''] )

workflow test_famas_trim {
    
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    FAMAS_TRIM ( input )
}
