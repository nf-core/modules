#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQ_ENA_PREP_ADAPTERREMOVAL_MD5SUM } from '../../../../subworkflows/nf-core/fastq_ena_prep_adapterremoval_md5sum/main.nf'


    input = [ [
                [ id:'test_se', single_end:true, collapse:false ], // meta map
                [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ],
                33
            ],
            [
                [ id:'test_pe', single_end:false, collapse:false ], // meta map
                [
                    file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                    file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
                ],
                33
            ] ]


workflow test_fastq_ena_prep_adapterremoval_md5sum {

    ch_input = Channel.fromList( input )

    FASTQ_ENA_PREP_ADAPTERREMOVAL_MD5SUM ( ch_input, [] )
}
