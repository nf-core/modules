#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_SORTVCF } from '../../../../../modules/nf-core/picard/sortvcf/main.nf'

workflow test_picard_sortvcf {

    input = [ [ id:'test' ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
            ]

    fasta  = [ [ id:'genome' ],
                file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
            ]

    dict   = [ [ id:'genome' ],
                file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)
            ]

    PICARD_SORTVCF ( input, fasta, dict )
}
