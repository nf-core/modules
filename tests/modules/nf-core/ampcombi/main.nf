#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AMPCOMBI } from '../../../../modules/nf-core/ampcombi/main.nf'
include { UNTAR as UNTAR1 ; UNTAR as UNTAR2 } from '../../../../modules/nf-core/untar/main.nf'

workflow test_ampcombi_directory {
   amp_input = [
       [ id:'sample_1'],
       file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files.tar.gz", checkIfExists: true)
   ]
   faa_input = [
       [ id:[]],
       file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_faa.tar.gz", checkIfExists: true)
   ]
   db = params.db ? file( params.db ) : []

   UNTAR1 ( amp_input )
   UNTAR2 ( faa_input )
   AMPCOMBI ( UNTAR1.out.untar, UNTAR2.out.untar.map{ it[1] }, db )
}

workflow test_ampcombi_file_paths {
    amp_input = [[id:'sample_1'],
                 [file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/ampir/sample_1/sample_1.ampir.tsv', checkIfExists: true),
                  file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/amplify/sample_1/sample_1_amplify.tsv', checkIfExists: true)]
                  ]
    faa_input = file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_faa/sample_1.faa', checkIfExists: true)
    
    db = params.db ? file( params.db ) : []

    AMPCOMBI ( amp_input, faa_input, db )
}