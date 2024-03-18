#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SVTK_COUNTSVTYPES } from '../../../../../modules/nf-core/svtk/countsvtypes/main.nf'

workflow test_svtk_countsvtypes {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true)
    ]

    SVTK_COUNTSVTYPES ( input )
}
