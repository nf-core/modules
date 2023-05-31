#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IPHOP_DOWNLOAD    } from '../../../../../modules/nf-core/iphop/download/main.nf'
include { IPHOP_PREDICT     } from '../../../../../modules/nf-core/iphop/predict/main.nf'

workflow test_iphop_predict {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://bitbucket.org/srouxjgi/iphop/raw/d27b6bbdcd39a6a1cb8407c44ccbcc800d2b4f78/test/test_input_phages.fna")
    ]

    IPHOP_DOWNLOAD ( )

    IPHOP_PREDICT ( input, IPHOP_DOWNLOAD.out.iphop_db )

}
