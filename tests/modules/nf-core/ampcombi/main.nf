#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AMPCOMBI } from '../../../../modules/nf-core/ampcombi/main.nf'

workflow test_ampcombi_file_paths {
    amp_input = [[id:'sample_1'],
                 [file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/ampir/sample_1/sample_1.ampir.tsv', checkIfExists: true),
                  file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/amplify/sample_1/sample_1_amplify.tsv', checkIfExists: true)]
                  ]
    faa_input = file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_faa/sample_1.faa', checkIfExists: true)
    
    db = params.db ? file( params.db ) : []

    AMPCOMBI ( amp_input, faa_input, db )
}