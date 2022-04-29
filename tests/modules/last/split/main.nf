#!/usr/bin/env nextflow



include { LAST_SPLIT } from '../../../../modules/last/split/main.nf'

workflow test_last_split {
    
    input = [ [ id:'contigs.genome' ], // meta map
              file(params.test_data['sarscov2']['genome']['contigs_genome_maf_gz'], checkIfExists: true) ]

    LAST_SPLIT ( input )
}
