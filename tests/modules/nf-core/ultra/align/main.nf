#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP        } from '../../../../../modules/nf-core/gunzip/main.nf'
include { GNU_SORT      } from '../../../../../modules/nf-core/gnu/sort/main.nf'
include { ULTRA_INDEX   } from '../../../../../modules/nf-core/ultra/index/main.nf'
include { ULTRA_ALIGN   } from '../../../../../modules/nf-core/ultra/align/main.nf'

workflow test_ultra_align {

    input_gunzip = [
        [ id:'test', single_end:false ],
        file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true)
    ]

    input_sort = [
        [ id:'test', single_end:false ],
        file(params.test_data['homo_sapiens']['genome']['genome_gtf']  , checkIfExists: true)
    ]

    genome = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    GUNZIP ( input_gunzip )
    GNU_SORT ( input_sort )
    ULTRA_INDEX ( genome, GNU_SORT.out.sorted.map{ it[1] } )
    ULTRA_ALIGN ( GUNZIP.out.gunzip, genome, ULTRA_INDEX.out.index )
}
