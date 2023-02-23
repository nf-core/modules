#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_FILTERINTERVALS } from '../../../../../modules/nf-core/gatk4/filterintervals/main.nf'

workflow test_gatk4_filterintervals {
    
    intervals = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_preprocessed_interval_list'], checkIfExists: true)
    ]
    read_counts = file(params.test_data['homo_sapiens']['genome']['genome_preprocessed_count_tsv'], checkIfExists: true)
    annotated_intervals = file(params.test_data['homo_sapiens']['genome']['genome_annotated_interval_tsv'], checkIfExists: true)

    GATK4_FILTERINTERVALS ( intervals, read_counts, annotated_intervals )
}
