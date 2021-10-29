#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IDR } from '../../../modules/idr/main.nf' addParams( options: [:] )

myFile = file('http://jordan.biology.gatech.edu/page/software/broadpeak/downloads/H3K27me3.bed', checkIfExists: true)
myFile.copyTo('peak_test')

workflow test_idr_narrowpeak {

    input = [
        file(params.test_data['homo_sapiens']['illumina']['test_narrow_peaks'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_narrow_peaks'], checkIfExists: true)
    ]

    IDR ( input, 'narrowPeak', 'test' )
}

workflow test_idr_broadpeak {

    input = [
        file(params.test_data['homo_sapiens']['illumina']['test_broad_peaks'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_broad_peaks'], checkIfExists: true)
    ]

    IDR ( input, 'broadPeak', 'test' )
}

workflow test_idr_noprefix {

    input = [
        file(params.test_data['homo_sapiens']['illumina']['test_narrow_peaks'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_narrow_peaks'], checkIfExists: true)
    ]

    IDR ( input, 'narrowPeak', '' )
}
