#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ALLELECOUNTER } from '../../../software/allelecounter/main.nf' addParams( options: [:] )

workflow test_allelecounter {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]
    positions = [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) ]

    ALLELECOUNTER ( input, positions )
}
