#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UCSC_WIGTOBIGWIG } from '../../../../software/ucsc/wigtobigwig/main.nf' addParams( options: [:] )

workflow test_ucsc_wigtobigwig {

    input = file(params.test_data['sarscov2']['illumina']['test_wig_gz'], checkIfExists: true)

    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)

    UCSC_WIGTOBIGWIG ( input, sizes )
}
