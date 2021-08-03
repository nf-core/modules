#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FUSIONCATCHER } from '../../../../modules/fusioncatcher/fusioncatcher/main.nf' addParams( options: [:] )

workflow test_fusioncatcher {

    input = [ [ id:'test', single_end:true ], // meta map
                [   file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_1_fastq_gz'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_2_fastq_gz'], checkIfExists: true) ]
            ]

    data_dir = file("/apps/ref/human_v102")
    FUSIONCATCHER ( input, data_dir )
}
