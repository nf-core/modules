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

// workflow test_idr_bed {

//     // Also tried with these two replicates in bed format from this course
//     // (http://jvanheld.github.io/cisreg_course/chip-seq/practical/annotation.html) and did not work
//     // http://pedagogix-tagc.univ-mrs.fr/courses/data/ngs/td_chip_seq/all/siGATA_ER_E2_r1_SRX176857_peaks.bed
//     // http://pedagogix-tagc.univ-mrs.fr/courses/data/ngs/td_chip_seq/all/siGATA_ER_E2_r2_SRX176859_peaks.bed

//     input = [
//         //Does not work with the data in test-datasets
//         // file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
//         // file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
//         // Also tried to see weather it swallow the same file as narrow and broad mode, but it does not
//         // file('https://raw.githubusercontent.com/kundajelab/idr/master/tests/data/peak1', checkIfExists: true),
//         // file('https://raw.githubusercontent.com/kundajelab/idr/master/tests/data/peak2', checkIfExists: true)
//     ]

//     IDR ( input, 'bed', 'test' )
// }
