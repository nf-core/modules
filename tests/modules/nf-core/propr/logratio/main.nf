#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PROPR_LOGRATIO } from '../../../../../modules/nf-core/propr/logratio/main.nf'


// workflow test_propr_logratio_default {

//     input = [
//         [ id:'test' ], // meta map
//         file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
//     ]

//     PROPR_LOGRATIO ( input )
// }

// workflow test_propr_logratio_clr_boxcox {

//     input = [
//         [ id:'test' ], // meta map
//         file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
//     ]

//     PROPR_LOGRATIO ( input )
// }

// workflow test_propr_logratio_alr {
    
//     input = [
//         [ id:'test' ], // meta map
//         file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
//     ]

//     PROPR_LOGRATIO ( input )
// }

// workflow test_propr_logratio_alr_boxcox {
    
//     input = [
//         [ id:'test' ], // meta map
//         file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
//     ]

//     PROPR_LOGRATIO ( input )
// }

// workflow test_propr_logratio_alr_geneid {
    
//     input = [
//         [ id:'test' ], // meta map
//         file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
//     ]

//     PROPR_LOGRATIO ( input )
// }

// workflow test_propr_logratio_alr_genename {
    
//     input = [
//         [ id:'test' ], // meta map
//         file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
//     ]

//     PROPR_LOGRATIO ( input )
// }

workflow test_propr_logratio_alr_geneid_boxcox {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    ]

    PROPR_LOGRATIO ( input )
}
