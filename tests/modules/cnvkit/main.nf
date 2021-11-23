#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT } from '../../../modules/cnvkit/main.nf'

workflow test_cnvkit {
    tumourbam = file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    normalbam = file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true)
    
    input = [ [ id:'test' ], // meta map
              tumourbam, 
              normalbam 
            ]
    fasta   = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    targets = file(params.test_data['sarscov2']['genome']['baits_bed'], checkIfExists: true)

    CNVKIT ( input, fasta, targets )
}
