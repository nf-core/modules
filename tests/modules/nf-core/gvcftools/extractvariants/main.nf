#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GVCFTOOLS_EXTRACTVARIANTS } from '../../../../../modules/nf-core/gvcftools/extractvariants/main.nf'

workflow test_gvcftools_extractvariants {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true)
    ]

    GVCFTOOLS_EXTRACTVARIANTS ( input )
}

workflow test_gvcftools_extractvariants_uncompressed {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)
    ]

    GVCFTOOLS_EXTRACTVARIANTS ( input )
}
