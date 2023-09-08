#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR } from '../../../../modules/nf-core/untar/main.nf'

workflow test_untar {
    input = [
        [],
        file(params.test_data['sarscov2']['genome']['kraken2_tar_gz'], checkIfExists: true)
    ]

    UNTAR ( input )
}

workflow test_untar_different_output_path {
    input = [
        [],
        file(params.test_data['homo_sapiens']['illumina']['test_flowcell'], checkIfExists: true)
    ]

    UNTAR ( input )
}


workflow test_untar_onlyfiles {
    input = [
        [],
        file(params.test_data['generic']['tar']['tar_gz'], checkIfExists: true)
    ]

    UNTAR ( input )
}
