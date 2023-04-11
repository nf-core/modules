#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ART_ILLUMINA } from '../../../../../modules/nf-core/art/illumina/main.nf'

workflow test_art_illumina_single_end {
    
    input = [ [ id: 'test'], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) 
    ]

    sequencing_system = 'HS25'
    fold_coverage = '15'
    read_length = '150'

    ART_ILLUMINA ( input, sequencing_system, fold_coverage, read_length  )
}

workflow test_art_illumina_paired_end {
    
    input = [ [ id: 'test'], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) 
    ]

    sequencing_system = 'HS25'
    fold_coverage = '15'
    read_length = '150'

    ART_ILLUMINA ( input, sequencing_system, fold_coverage, read_length  )
}