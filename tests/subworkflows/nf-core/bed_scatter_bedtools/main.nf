#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BED_SCATTER_BEDTOOLS } from '../../../../subworkflows/nf-core/bed_scatter_bedtools/main.nf'

workflow test_bed_scatter_bedtools_bed {

    input = Channel.of([
        [id:'test'],
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true),
        2
    ])

    BED_SCATTER_BEDTOOLS (
        input
    )
}
