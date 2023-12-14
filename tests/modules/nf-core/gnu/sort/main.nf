#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { GNU_SORT as GNU_SORT_A } from '../../../../../modules/nf-core/gnu/sort/main.nf'
include { GNU_SORT as GNU_SORT_B } from '../../../../../modules/nf-core/gnu/sort/main.nf'
include { GNU_SORT as GNU_SORT_C } from '../../../../../modules/nf-core/gnu/sort/main.nf'

workflow test_gnu_sort {
    
    input_a = [
        [id:'test'],
        file(params.test_data['generic']['unsorted_data']['unsorted_text']['genome_file'],
        checkIfExists: true)
    ]

    GNU_SORT_A ( input_a )

    input_b = [
        [id:'test'],
        file(params.test_data['generic']['unsorted_data']['unsorted_text']['intervals'],
        checkIfExists: true)
    ]

    GNU_SORT_B ( input_b )

    input_c = [
        [id:'test'],
        file(params.test_data['generic']['unsorted_data']['unsorted_text']['numbers_csv'],
        checkIfExists: true)
    ]

    GNU_SORT_C ( input_c )

}
