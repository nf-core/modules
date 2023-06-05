#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_MARKDUPLICATES as PICARD_MARKDUPLICATES_SORTED_BAM } from '../../../../../modules/nf-core/picard/markduplicates/main.nf'
include { PICARD_MARKDUPLICATES as PICARD_MARKDUPLICATES_UNSORTED_BAM}  from '../../../../../modules/nf-core/picard/markduplicates/main.nf'
include { PICARD_MARKDUPLICATES as PICARD_MARKDUPLICATES_SORTED_CRAM }  from '../../../../../modules/nf-core/picard/markduplicates/main.nf'

workflow test_picard_markduplicates_sorted_bam  {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]
    fasta = [ [ id:'genome' ],
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
            ]
    fai   = [ [ id:'genome' ],
              file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
            ]

    PICARD_MARKDUPLICATES_SORTED_BAM ( input, fasta, fai )
}

workflow test_picard_markduplicates_unsorted_bam  {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
            ]
    fasta = [ [ id:'genome' ],
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
            ]
    fai   = [ [ id:'genome' ],
              file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
            ]

    PICARD_MARKDUPLICATES_UNSORTED_BAM ( input, fasta, fai )
}

workflow test_picard_markduplicates_sorted_cram  {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true)
            ]
    fasta = [ [ id:'genome' ],
              file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
            ]
    fai   = [ [ id:'genome' ],
              file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
            ]

    PICARD_MARKDUPLICATES_SORTED_CRAM ( input, fasta, fai )
}

// TODO Add a test with PICARD_MARKDUPLICATES ( input, [], [] )
