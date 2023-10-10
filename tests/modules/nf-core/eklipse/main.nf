#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EKLIPSE } from '../../../../modules/nf-core/eklipse/main.nf'

workflow test_eklipse {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_illumina_mt_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_illumina_mt_bam_bai'], checkIfExists: true)
    ]
    ref_gb = [ file(params.test_data['homo_sapiens']['genome']['genome_mt_gb'], checkIfExists: true) ]

    EKLIPSE ( input, ref_gb )
}
