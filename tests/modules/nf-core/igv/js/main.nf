#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IGV_JS } from '../../../../../modules/nf-core/igv/js/main.nf'

workflow test_igv_js {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    IGV_JS ( input )
}
