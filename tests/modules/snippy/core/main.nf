#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNIPPY_RUN } from '../../../../modules/snippy/run/main.nf'
include { SNIPPY_CORE } from '../../../../modules/snippy/core/main.nf'

workflow test_snippy_core {
    
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['candidatus_portiera_aleyrodidarum']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['candidatus_portiera_aleyrodidarum']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    reference = file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['genome_fasta'], checkIfExists: true)

    SNIPPY_RUN ( input, reference )
    SNIPPY_RUN.out.vcf.join(SNIPPY_RUN.out.aligned_fa).map{ meta, vcf, aln -> [[id:'snippy-core'], vcf, aln] }.set{ ch_snippy_core }
    SNIPPY_CORE(ch_snippy_core, reference)

}
