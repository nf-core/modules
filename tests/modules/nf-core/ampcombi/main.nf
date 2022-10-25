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

//workflow test_ampcombi_file_paths {
//    amp_input = 
//        channel.fromPath(
//       [ 
//        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/ampir/sample_1/sample_1.ampir.tsv'),
//        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/amplify/sample_1/sample_1_amplify.tsv')
//        ] , checkIfExists: true )
//        .collect { it as String }
//        .map { 
//            def fmeta = ["id": "sample_1"]
//            [ fmeta, it ]
//            }
//    //    [ id:'sample_1' ],
//    //    ['https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/ampir/sample_1/sample_1.ampir.tsv',
//    //    'https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/amplify/sample_1/sample_1_amplify.tsv']
//    //]
//    faa_folder = [
//        [ id:'test' ],
//        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_faa.tar.gz", checkIfExists: true)
//    ]
//    
//    UNTAR2 ( faa_folder )
//    AMPCOMBI ( amp_input, UNTAR2.out.untar.map{ it[1] } )
//}
workflow test_ampcombi_file_paths {
    amp_input = [[id:'sample_1'],
                 [file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/ampir/sample_1/sample_1.ampir.tsv', checkIfExists: true),
                  file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/amplify/sample_1/sample_1_amplify.tsv', checkIfExists: true)]]

    faa_folder = [
        [ id:'sample_1' ],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_faa.tar.gz", checkIfExists: true)
    ]
    
    UNTAR2 ( faa_folder )
    AMPCOMBI ( amp_input, UNTAR2.out.untar.map{ it[1] } )
}


//workflow test_ampcombi_file_paths {
//    ampir_input_ch = Channel.from(
//            [ id:'sample_1' ],
//            file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/ampir/sample_1/sample_1.ampir.tsv')
//        )
//    amplify_input_ch = Channel.from(
//           [ id:'sample_1' ],
//            file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/amplify/sample_1/sample_1_amplify.tsv')
//        )
//    
//    amp_input = 

    //    [ id:'sample_1' ],
    //    ['https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/ampir/sample_1/sample_1.ampir.tsv',
    //    'https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_files/amplify/sample_1/sample_1_amplify.tsv']
    //]
//    faa_folder = [
//        [ id:'test' ],
 //       file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/ampcombi/test_faa.tar.gz", checkIfExists: true)
 //   ]
 //   
 //   UNTAR2 ( faa_folder )
 //   AMPCOMBI ( amp_input, UNTAR2.out.untar.map{ it[1] } )
//}