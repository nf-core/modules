#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNIPPY_RUN } from '../../../../../modules/nf-core/snippy/run/main.nf'
include { SNIPPY_CORE } from '../../../../../modules/nf-core/snippy/core/main.nf'

workflow test_snippy_core {
    
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['candidatus_portiera_aleyrodidarum']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['candidatus_portiera_aleyrodidarum']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    reference = file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['genome_fasta'], checkIfExists: true)

    SNIPPY_RUN ( input, reference )
    SNIPPY_RUN.out.vcf.collect{meta, vcf -> vcf}.map{ vcf -> [[id:'snippy-core'], vcf]}.set{ ch_merge_vcf }
    SNIPPY_RUN.out.aligned_fa.collect{meta, aligned_fa -> aligned_fa}.map{ aligned_fa -> [[id:'snippy-core'], aligned_fa]}.set{ ch_merge_aligned_fa }
    ch_merge_vcf.join( ch_merge_aligned_fa ).set{ ch_snippy_core }
    SNIPPY_CORE( ch_snippy_core, reference )

}
