#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_REHEADER } from '../../../../modules/bcftools/reheader/main.nf' addParams( options: [suffix: '.updated'] )

workflow test_bcftools_reheader_sequences {

    input = [
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true)
    ]
    fai    = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    header = []


    BCFTOOLS_REHEADER ( input, fai, [] )
}
