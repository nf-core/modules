#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_FIXREF } from '../../../../modules/bcftools/fixref/main.nf' addParams( options: [:] )

workflow test_bcftools_fixref {
    
    input = [ [ id:'test2', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true)]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)


    BCFTOOLS_FIXREF ( input )
}
