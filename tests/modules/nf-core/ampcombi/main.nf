#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AMPCOMBI } from '../../../../modules/nf-core/ampcombi/main.nf'
include { UNTAR as UNTAR1 ; UNTAR as UNTAR2 } from '../../../../modules/nf-core/untar/main.nf'

workflow test_ampcombi {
    amp_input = [
        [ id:'test' ],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files.tar.gz", checkIfExists: true)
    ]
    faa_folder = [
        [ id:'test' ],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_faa.tar.gz", checkIfExists: true)
    ]
    outdir = "ampcombi_results"

    UNTAR1 ( amp_input )
    UNTAR2 ( faa_folder )
    AMPCOMBI ( UNTAR1.out.untar, UNTAR2.out.untar.map{ it[1] }, outdir )
}