#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ASHLAR } from '../../../../modules/nf-core/ashlar/main.nf'
include { ASHLAR as ASHLAR_TILE } from '../../../../modules/nf-core/ashlar/main.nf'

workflow test_ashlar_1_file {

    input_list =  [ [ id:'test_all' ],
               [file(params.test_data['imaging']
                                     ['ome-tiff']
                                     ['cycif_tonsil_cycle1'], checkIfExists: true)] ]

    ASHLAR ( input_list, [], [] )

}

workflow test_ashlar_all_files {

    input_list =  [ [ id:'test_all' ],
               [file(params.test_data['imaging']
                                     ['ome-tiff']
                                     ['cycif_tonsil_cycle1'], checkIfExists: true),
                file(params.test_data['imaging']
                                     ['ome-tiff']
                                     ['cycif_tonsil_cycle2'], checkIfExists: true),
                file(params.test_data['imaging']
                                     ['ome-tiff']
                                     ['cycif_tonsil_cycle3'], checkIfExists: true)] ]

    ASHLAR ( input_list, [], [] )

}

workflow test_ashlar_all_files_tile_size {

    input_list =  [ [ id:'test_all' ],
               [file(params.test_data['imaging']
                                     ['ome-tiff']
                                     ['cycif_tonsil_cycle1'], checkIfExists: true),
                file(params.test_data['imaging']
                                     ['ome-tiff']
                                     ['cycif_tonsil_cycle2'], checkIfExists: true),
                file(params.test_data['imaging']
                                     ['ome-tiff']
                                     ['cycif_tonsil_cycle3'], checkIfExists: true)] ]

    ASHLAR_TILE ( input_list, [], [] )

}

workflow test_ashlar_all_files_dfp_ffp {

    input_list =  [ [ id:'test_all' ],
               [file(params.test_data['imaging']
                                     ['ome-tiff']
                                     ['cycif_tonsil_cycle1'], checkIfExists: true),
                file(params.test_data['imaging']
                                     ['ome-tiff']
                                     ['cycif_tonsil_cycle2'], checkIfExists: true),
                file(params.test_data['imaging']
                                     ['ome-tiff']
                                     ['cycif_tonsil_cycle3'], checkIfExists: true)] ]

    ASHLAR ( input_list, [file(params.test_data['imaging']
                                               ['ome-tiff']
                                               ['cycif_tonsil_dfp'], checkIfExists: true)],
                         [file(params.test_data['imaging']
                                               ['ome-tiff']
                                               ['cycif_tonsil_ffp'], checkIfExists: true)] )

}

workflow test_ashlar_1_file_dfp_ffp {

    input_list =  [ [ id:'test_all' ],
               [file(params.test_data['imaging']
                                     ['ome-tiff']
                                     ['cycif_tonsil_cycle1'], checkIfExists: true)] ]

    ASHLAR ( input_list, [file(params.test_data['imaging']
                                               ['ome-tiff']
                                               ['cycif_tonsil_dfp'], checkIfExists: true)],
                         [file(params.test_data['imaging']
                                               ['ome-tiff']
                                               ['cycif_tonsil_ffp'], checkIfExists: true)] )

}

