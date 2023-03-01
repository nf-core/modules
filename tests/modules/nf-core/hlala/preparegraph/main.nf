#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HLALA_PREPAREGRAPH } from '../../../../../modules/nf-core/hlala/preparegraph/main.nf'
include { UNZIP } from '../../../../../modules/nf-core/unzip/main.nf'

workflow test_hlala_preparegraph {
    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['genome']['prg_input'], checkIfExists: true)
    ]
    
    UNZIP( input )

    UNZIP.out.unzipped_archive.dump(tag: "unzipped archive")

    HLALA_PREPAREGRAPH ( UNZIP.out.unzipped_archive )
}
