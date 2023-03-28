#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR_FILES } from '../../../../modules/nf-core/untar_files/main.nf'

workflow test_untar_files {
    input = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['kraken2_tar_gz'], checkIfExists: true)
    ]

    UNTAR_FILES ( input )
}

workflow test_untar_files_different_output_path {
    input = [
        [id: 'test'],
        file(params.test_data['homo_sapiens']['illumina']['test_flowcell'], checkIfExists: true)
    ]

    UNTAR_FILES ( input )
}


workflow test_untar_files_onlyfiles {
    input = [
        [id: 'test'],
        file(params.test_data['generic']['tar']['tar_gz'], checkIfExists: true)
    ]

    UNTAR_FILES ( input )
}
