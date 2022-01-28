#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UCSC_WIGTOBIGWIG } from '../../../../modules/ucsc/wigtobigwig/main.nf'

workflow test_ucsc_wigtobigwig {

    input = [ [ id:'test', single_end:false ], // meta map,
                file(params.test_data['sarscov2']['illumina']['test_wig_gz'], checkIfExists: true) ]

    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)

    UCSC_WIGTOBIGWIG ( input, sizes )
}
