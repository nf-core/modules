#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_BAMTOBED } from '../../../../software/bedtools/bamtobed/main.nf' addParams( options: [:] )

workflow test_bedtools_bamtobed {
    input = [ [ id:'test'], //meta map
              file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true),
              file(params.test_data['sarscov2']['nanopore']['test_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['nanopore']['test_sorted_bam_bai'], checkIfExists: true)
            ]

    BEDTOOLS_BAMTOBED ( input )
}
