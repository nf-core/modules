#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AMPCOMBI } from '../../../../modules/nf-core/ampcombi/main.nf'
include { UNTAR as UNTAR1 ; UNTAR as UNTAR2 } from '../../../../modules/nf-core/untar/main.nf'

// workflow test_ampcombi_directory {
//    amp_input = [
//        [ id:'test'],
//        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files.tar.gz", checkIfExists: true)
//    ]
//    faa_folder = [
//        [ id:'test'],
//        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_faa.tar.gz", checkIfExists: true)
//    ]

//    UNTAR1 ( amp_input )
//    UNTAR2 ( faa_folder )
//    AMPCOMBI ( UNTAR1.out.untar, UNTAR2.out.untar.map{ it[1] } )
// }

workflow test_ampcombi_file_paths {
    amp_input = [
    //    channel.fromPath(
    //   [ 
    //        'https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/ampir/sample_1/sample_1.ampir.tsv',
    //        'https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/amplify/sample_1/sample_1_amplify.tsv',
    //    ] , checkIfExists: true ).collect()
        [ id:'sample_1' ],
        ['https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/ampir/sample_1/sample_1.ampir.tsv',
        'https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/amplify/sample_1/sample_1_amplify.tsv']
    ]

    faa_folder = [
        [ id:'test' ],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_faa.tar.gz", checkIfExists: true)
    ]
    
    UNTAR2 ( faa_folder )
    AMPCOMBI ( amp_input, UNTAR2.out.untar.map{ it[1] } )
}