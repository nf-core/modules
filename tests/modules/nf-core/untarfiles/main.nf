#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTARFILES } from '../../../../modules/nf-core/untarfiles/main.nf'

workflow test_untarfiles {
    input = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['kraken2_tar_gz'], checkIfExists: true)
    ]

    UNTARFILES ( input )
}

workflow test_untarfiles_different_output_path {
    input = [
        [id: 'test'],
        file(params.test_data['homo_sapiens']['illumina']['test_flowcell'], checkIfExists: true)
    ]

    UNTARFILES ( input )
}


workflow test_untarfiles_onlyfiles {
    input = [
        [id: 'test'],
        file(params.test_data['generic']['tar']['tar_gz'], checkIfExists: true)
    ]

    UNTARFILES ( input )
}
