#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METAEUK_EASYPREDICT } from '../../../../../modules/nf-core/metaeuk/easypredict/main.nf'

workflow test_metaeuk_easypredict {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    database = [
        file(params.test_data['proteomics']['database']['yeast_ups'], checkIfExists: true)
    ]

    METAEUK_EASYPREDICT ( input, database )
}
