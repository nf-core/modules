#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IDR } from '../../../modules/idr/main.nf'

workflow test_idr_narrowpeak {

    input = [
        file(params.test_data['homo_sapiens']['illumina']['test_narrowpeak'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_narrowpeak'], checkIfExists: true)
    ]

    IDR ( input, 'narrowPeak', 'test' )
}

workflow test_idr_broadpeak {

    input = [
        file(params.test_data['homo_sapiens']['illumina']['test_broadpeak'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_broadpeak'], checkIfExists: true)
    ]

    IDR ( input, 'broadPeak', 'test' )
}

workflow test_idr_noprefix {

    input = [
        file(params.test_data['homo_sapiens']['illumina']['test_narrowpeak'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_narrowpeak'], checkIfExists: true)
    ]

    IDR ( input, 'narrowPeak', '' )
}
