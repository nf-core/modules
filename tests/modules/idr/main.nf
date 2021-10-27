#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IDR } from '../../../modules/idr/main.nf' addParams( options: [:] )

// TODO: Add appropriate test data to nf-core/testdatasets and change paths below

workflow test_idr_narrowpeak {

    input = [
        //file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
        //file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
        file('https://raw.githubusercontent.com/kundajelab/idr/master/tests/data/peak1', checkIfExists: true),
        file('https://raw.githubusercontent.com/kundajelab/idr/master/tests/data/peak2', checkIfExists: true)
    ]

    IDR ( input, 'narrowPeak', 'test' )
}

// workflow test_idr_broadpeak {

//     input = [
//         //file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
//         //file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
//         file('https://raw.githubusercontent.com/kundajelab/idr/master/tests/data/peak1', checkIfExists: true),
//         file('https://raw.githubusercontent.com/kundajelab/idr/master/tests/data/peak2', checkIfExists: true)
//     ]

//     IDR ( input, 'broadPeak', 'test' )
// }

// workflow test_idr_bed {

//     input = [
//         //file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
//         //file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
//         file('https://raw.githubusercontent.com/kundajelab/idr/master/tests/data/peak1', checkIfExists: true),
//         file('https://raw.githubusercontent.com/kundajelab/idr/master/tests/data/peak2', checkIfExists: true)
//     ]

//     IDR ( input, 'bed', 'test' )
// }

workflow test_idr_noprefix {

    input = [
        //file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
        //file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
        file('https://raw.githubusercontent.com/kundajelab/idr/master/tests/data/peak1', checkIfExists: true),
        file('https://raw.githubusercontent.com/kundajelab/idr/master/tests/data/peak2', checkIfExists: true)
    ]

    IDR ( input, 'narrowPeak', '' )
}
