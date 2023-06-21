#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ASHLAR } from '../../../../modules/nf-core/ashlar/main.nf'
// we zero out the UUID of output tiff images with ZERO_UUID so we get a consistent md5sum
include { ZERO_UUID } from './zero_uuid.nf'

workflow test_ashlar_1_file {

    input_list =  [ [ [ id:'test_all' ],
               [file(params.test_data['imaging']['ome-tiff']['cycif_tonsil_cycle1'], checkIfExists: true)] ] ]
    input_channel = Channel.fromList(input_list)

    ASHLAR ( input_channel )

    ZERO_UUID ( ASHLAR.out[0], "8390123" )

}

workflow test_ashlar_all_files {

    input_list =  [ [ [ id:'test_all' ],
               [file(params.test_data['imaging']['ome-tiff']['cycif_tonsil_cycle1'], checkIfExists: true),
                file(params.test_data['imaging']['ome-tiff']['cycif_tonsil_cycle2'], checkIfExists: true),
                file(params.test_data['imaging']['ome-tiff']['cycif_tonsil_cycle3'], checkIfExists: true)] ] ]
    input_channel = Channel.fromList(input_list)

    ASHLAR ( input_channel )

    ZERO_UUID ( ASHLAR.out[0], "25169643" )

}
