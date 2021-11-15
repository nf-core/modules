#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQSCAN } from '../../../modules/fastqscan/main.nf' addParams( options: [ args: "-g 30000"] )

workflow test_fastqscan {
    
    input = [ [ id:'test', single_end:true ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]

    FASTQSCAN ( input )
}
