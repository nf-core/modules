#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQ_ENA_PREP_ADAPTERREMOVAL } from '../../../../subworkflows/nf-core/fastq_ena_prep_adapterremoval/main.nf'


    input = [ [
                [ id:'test', single_end:true, collapse:false ], // meta map
                [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ],
                33
            ],
            [
                [ id:'test', single_end:false, collapse:false ], // meta map
                [
                    file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                    file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
                ],
                33
            ] ]


workflow test_fastq_ena_prep_adapterremoval {

    ch_input = Channel.fromList( input )
        .dump(tag:'input')

    FASTQ_ENA_PREP_ADAPTERREMOVAL ( ch_input, [] )
}
