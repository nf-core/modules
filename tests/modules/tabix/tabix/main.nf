#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TABIX_TABIX as TABIX_BED } from '../../../../modules/tabix/tabix/main.nf'
include { TABIX_TABIX as TABIX_GFF } from '../../../../modules/tabix/tabix/main.nf'
include { TABIX_TABIX as TABIX_VCF } from '../../../../modules/tabix/tabix/main.nf'

workflow test_tabix_tabix_bed {
    input = [ [ id:'B.bed' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['test_bed_gz'], checkIfExists: true) ] 
            ]

    TABIX_BED ( input )
}

workflow test_tabix_tabix_gff {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['genome_gff3_gz'], checkIfExists: true) ] 
            ]

    TABIX_GFF ( input )
}

workflow test_tabix_tabix_vcf {
    input = [ [ id:'test.vcf' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true) ] 
            ]

    TABIX_VCF ( input )
}
