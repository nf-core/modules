#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

test_options = ['args': '--filter-name "test_filter" --filter-expression "MQ0 > 0"', 'suffix': '.filtered']
include { GATK4_VARIANTFILTRATION } from '../../../../software/gatk4/variantfiltration/main.nf' addParams( options: test_options )

workflow test_gatk4_variantfiltration {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]
    fasta = [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
    fai = [ file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true) ]
    genome_dict = [ file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true) ]


    GATK4_VARIANTFILTRATION ( input, fasta, fai, genome_dict )
}
