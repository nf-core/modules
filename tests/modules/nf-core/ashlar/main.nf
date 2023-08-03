#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ASHLAR } from '../../../../modules/nf-core/ashlar/main.nf'
include { ASHLAR as ASHLAR_TILE } from '../../../../modules/nf-core/ashlar/main.nf'

// we zero out the UUID of output tiff images with ZERO_UUID so we get a consistent md5sum
process ZERO_UUID {

    input:
    val(file_in)
    val(offset)

    when:
    file_in != "versions.yml"

    script:
    def file_path = file_in[1]

    """
    echo -n "00000000-0000-0000-0000-000000000000" | dd of=$file_path bs=1 seek=$offset conv=notrunc
    """

}

workflow test_ashlar_1_file {

    input_list =  [ [ id:'test_all' ],
               [file(params.test_data['imaging']
                                     ['ome-tiff']
                                     ['cycif_tonsil_cycle1'], checkIfExists: true)] ]

    ASHLAR ( input_list, [], [] )

    ZERO_UUID ( ASHLAR.out.tif, "8390123" )

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

    ZERO_UUID ( ASHLAR.out.tif, "25169643" )

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

    ZERO_UUID ( ASHLAR_TILE.out.tif, "12586923" )

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

    ZERO_UUID ( ASHLAR.out.tif, "25169643" )

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

    ZERO_UUID ( ASHLAR.out.tif, "8390123" )

}

