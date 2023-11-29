#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TABIX_TABIX as TABIX_BED     } from '../../../../../modules/nf-core/tabix/tabix/main.nf'
include { TABIX_TABIX as TABIX_GFF     } from '../../../../../modules/nf-core/tabix/tabix/main.nf'
include { TABIX_TABIX as TABIX_VCF_TBI } from '../../../../../modules/nf-core/tabix/tabix/main.nf'
include { TABIX_TABIX as TABIX_VCF_CSI } from '../../../../../modules/nf-core/tabix/tabix/main.nf'

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

workflow test_tabix_tabix_vcf_tbi {
    input = [ [ id:'test.vcf' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true) ] 
            ]

    TABIX_VCF_TBI ( input )
}

workflow test_tabix_tabix_vcf_csi {
    input = [ [ id:'test.vcf' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true) ]
            ]

    TABIX_VCF_CSI ( input )
}
