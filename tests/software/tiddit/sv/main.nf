#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TIDDIT_SV } from '../../../../software/tiddit/sv/main.nf' addParams( options: [:] )

workflow test_tiddit_sv {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ] 
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    TIDDIT_SV ( input, fasta, fai )
}

workflow test_tiddit_sv_no_ref {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ] 
            ]

    TIDDIT_SV ( input, [], [] )
}
