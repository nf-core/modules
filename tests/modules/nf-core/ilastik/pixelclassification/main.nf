#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ILASTIK_PIXELCLASSIFICATION as ILASTIK_PIXELCLASSIFICATION_PROBABILITIES } from '../../../../../modules/nf-core/ilastik/pixelclassification/main.nf'
include { ILASTIK_PIXELCLASSIFICATION as ILASTIK_PIXELCLASSIFICATION_SIMPLESEGMENTATION } from '../../../../../modules/nf-core/ilastik/pixelclassification/main.nf'
include { ILASTIK_PIXELCLASSIFICATION as ILASTIK_PIXELCLASSIFICATION_UNCERTAINTY } from '../../../../../modules/nf-core/ilastik/pixelclassification/main.nf'
include { ILASTIK_PIXELCLASSIFICATION as ILASTIK_PIXELCLASSIFICATION_FEATURES } from '../../../../../modules/nf-core/ilastik/pixelclassification/main.nf'

workflow test_ilastik_pixelclassification {

    input = [
        [ id:'probabilities' ], // meta map
        file(params.test_data['imaging']['h5']['plant_wga'], checkIfExists: true)
    ]

    ilp   = [
        [id:'project'],
        file(params.test_data['imaging']['ilp']['plant_wga_pixel_class'], checkIfExists: true)
    ]

    ILASTIK_PIXELCLASSIFICATION_PROBABILITIES ( input, ilp )
}

workflow test_ilastik_pixelclassification_simplesegmentation {

    input = [
        [ id:'segmentations' ], // meta map
        file(params.test_data['imaging']['h5']['plant_wga'], checkIfExists: true)
    ]
    ilp   = [
        [id:'project'],
        file(params.test_data['imaging']['ilp']['plant_wga_pixel_class'], checkIfExists: true)
    ]

    ILASTIK_PIXELCLASSIFICATION_SIMPLESEGMENTATION ( input, ilp )
}

workflow test_ilastik_pixelclassification_uncertainty {

    input = [
        [ id:'uncertainty' ], // meta map
        file(params.test_data['imaging']['h5']['plant_wga'], checkIfExists: true)
    ]
    ilp   = [
        [id:'project'],
        file(params.test_data['imaging']['ilp']['plant_wga_pixel_class'], checkIfExists: true)
    ]

    ILASTIK_PIXELCLASSIFICATION_UNCERTAINTY ( input, ilp )
}

workflow test_ilastik_pixelclassification_features {

    input = [
        [ id:'features' ], // meta map
        file(params.test_data['imaging']['h5']['plant_wga'], checkIfExists: true)
    ]
    ilp   = [
        [id:'project'],
        file(params.test_data['imaging']['ilp']['plant_wga_pixel_class'], checkIfExists: true)
    ]

    ILASTIK_PIXELCLASSIFICATION_FEATURES ( input, ilp )
}
