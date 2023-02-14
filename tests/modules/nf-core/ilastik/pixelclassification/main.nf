#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ILASTIK_PIXELCLASSIFICATION as ILASTIK_PIXELCLASSIFICATION_PROBABILITIES } from '../../../../../modules/nf-core/ilastik/pixelclassification/main.nf'
include { ILASTIK_PIXELCLASSIFICATION as ILASTIK_PIXELCLASSIFICATION_SIMPLESEGMENTATION } from '../../../../../modules/nf-core/ilastik/pixelclassification/main.nf'
include { ILASTIK_PIXELCLASSIFICATION as ILASTIK_PIXELCLASSIFICATION_UNCERTAINTY } from '../../../../../modules/nf-core/ilastik/pixelclassification/main.nf'
include { ILASTIK_PIXELCLASSIFICATION as ILASTIK_PIXELCLASSIFICATION_FEATURES } from '../../../../../modules/nf-core/ilastik/pixelclassification/main.nf'

workflow test_ilastik_pixelclassification {

    input = [
        [ id:'test' ], // meta map
        file("/Users/florian_wuennemann/1_Projects/nf_core/test_data/wga.h5")
    ]
    ilp = [file("/Users/florian_wuennemann/1_Projects/nf_core/test_data/plant_wga.ilp")]

    ILASTIK_PIXELCLASSIFICATION_PROBABILITIES ( input, ilp )
}

workflow test_ilastik_pixelclassification_simplesegmentation {

    input = [
        [ id:'test' ], // meta map
        file("/Users/florian_wuennemann/1_Projects/nf_core/test_data/wga.h5")
    ]
    ilp = [file("/Users/florian_wuennemann/1_Projects/nf_core/test_data/plant_wga.ilp")]

    ILASTIK_PIXELCLASSIFICATION_SIMPLESEGMENTATION ( input, ilp )
}

workflow test_ilastik_pixelclassification_uncertainty {

    input = [
        [ id:'test' ], // meta map
        file("/Users/florian_wuennemann/1_Projects/nf_core/test_data/wga.h5")
    ]
    ilp = [file("/Users/florian_wuennemann/1_Projects/nf_core/test_data/plant_wga.ilp")]

    ILASTIK_PIXELCLASSIFICATION_UNCERTAINTY ( input, ilp )
}

workflow test_ilastik_pixelclassification_features {

    input = [
        [ id:'test' ], // meta map
        file("/Users/florian_wuennemann/1_Projects/nf_core/test_data/wga.h5")
    ]
    ilp = [file("/Users/florian_wuennemann/1_Projects/nf_core/test_data/plant_wga.ilp")]

    ILASTIK_PIXELCLASSIFICATION_FEATURES ( input, ilp )
}
