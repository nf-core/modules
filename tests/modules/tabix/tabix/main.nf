#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TABIX_TABIX as TABIX_BED } from '../../../../software/tabix/tabix/main.nf' addParams( options: ['args': '-p bed'] )
include { TABIX_TABIX as TABIX_GFF } from '../../../../software/tabix/tabix/main.nf' addParams( options: ['args': '-p gff'] )
include { TABIX_TABIX as TABIX_VCF } from '../../../../software/tabix/tabix/main.nf' addParams( options: ['args': '-p vcf'] )

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
