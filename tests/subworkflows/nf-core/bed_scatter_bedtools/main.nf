#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BED_SCATTER_BEDTOOLS } from '../../../../subworkflows/nf-core/bed_scatter_bedtools/main.nf'

workflow test_bed_scatter_bedtools_bed {

    input = Channel.of([
        [id:'test'],
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true),
        2
    ])

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    BED_SCATTER_BEDTOOLS (
        input,
        fasta_fai
    )
}

workflow test_bed_scatter_bedtools_no_bed {

    input = Channel.of([
        [id:'test'],
        [],
        2
    ])

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    BED_SCATTER_BEDTOOLS (
        input,
        fasta_fai
    )
}

workflow test_bed_scatter_bedtools_mixed {

    input = Channel.of([
        [id:'test'],
        [],
        2
    ],
    [
        [id:'test2'],
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true),
        2
    ])

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    BED_SCATTER_BEDTOOLS (
        input,
        fasta_fai
    )
}
