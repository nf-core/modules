#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FAMSA_GUIDETREE                           } from '../../../../../modules/nf-core/famsa/guidetree/main.nf'
include { TCOFFEE_ALIGN as TCOFFEE_ALIGN_SEQUENCE   } from '../../../../../modules/nf-core/tcoffee/align/main.nf'
include { TCOFFEE_ALIGN as TCOFFEE_ALIGN_WITHTREE   } from '../../../../../modules/nf-core/tcoffee/align/main.nf'
include { TCOFFEE_ALIGN as TCOFFEE_ALIGN_STRUCTURES } from '../../../../../modules/nf-core/tcoffee/align/main.nf'
include { UNTAR                                     } from '../../../../../modules/nf-core/untar/main.nf'


// workflow test_tcoffee_align_sequence {
    
//     input = [
//         [ id:'test' ],
//         file(params.test_data['sarscov2']['genome']['informative_sites_fas'], checkIfExists: true)
//     ]

//     TCOFFEE_ALIGN_SEQUENCE ( input,  [[:],[]],  [[:],[],[]] )
// }


// workflow test_famsa_align_with_tree {
    
//     input = [
//         [ id:'test' ],
//         file(params.test_data['sarscov2']['genome']['informative_sites_fas'], checkIfExists: true)
//     ]

    
//     ch_tree = FAMSA_GUIDETREE ( input ).tree

//     TCOFFEE_ALIGN_WITHTREE ( input , ch_tree,  [[:],[],[]])

// }





workflow test_famsa_align_with_structures {
    
    input = [
        [ id:'test' ],
        file("https://raw.githubusercontent.com/nf-core/test-datasets/multiplesequencealign/testdata/setoxin-ref.fa", checkIfExists: true)
    ]

    structures = [
        [ id:'test' ],
        file("https://raw.githubusercontent.com/nf-core/test-datasets/multiplesequencealign/testdata/structures/seatoxin-ref.tar.gz", checkIfExists: true)
    ]


    ch_structures = UNTAR ( structures ).untar
    ch_structures = ch_structures.map { meta,dir -> [[ id:'test' ], [] ,file(dir).listFiles().collect()]}

    ch_structures
            .map { meta,template,file ->
                def id = file.nameWithoutExtension
                "${id}_P_${id}"
            }
            .flatten()
            .splitText()
            .view()


    TCOFFEE_ALIGN_STRUCTURES ( input ,  [[:],[]] , ch_structures)

}