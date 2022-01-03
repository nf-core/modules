#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GSTAMA_MERGE } from '../../../../modules/gstama/merge/main'

workflow test_gstama_merge {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['pacbio']['genemodel1'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['pacbio']['genemodel2'], checkIfExists: true)
        ]
    ]
    filelist = file(params.test_data['homo_sapiens']['pacbio']['filelist'], checkIfExists: true)

    GSTAMA_MERGE ( input, filelist )
}
