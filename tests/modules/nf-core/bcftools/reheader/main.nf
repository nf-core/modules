#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_REHEADER } from '../../../../../modules/nf-core/bcftools/reheader/main.nf'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_VCF_GZ} from '../../../../../modules/nf-core/bcftools/reheader/main.nf'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_BCF_GZ} from '../../../../../modules/nf-core/bcftools/reheader/main.nf'

workflow test_bcftools_reheader_update_sequences {

    input = [
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        []
    ]
    fai    = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    BCFTOOLS_REHEADER ( input, fai )
}

workflow test_bcftools_reheader_update_sequences_vcf_gz {

    input = [
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        []
    ]
    fai    = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    BCFTOOLS_REHEADER_VCF_GZ ( input, fai )
}

workflow test_bcftools_reheader_update_sequences_bcf_gz {

    input = [
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        []
    ]
    fai    = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    BCFTOOLS_REHEADER_BCF_GZ ( input, fai )
}

workflow test_bcftools_reheader_new_header {

    input = [
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]
    fai    = []

    BCFTOOLS_REHEADER ( input, fai )
}

workflow test_bcftools_reheader_new_header_update_sequences {

    input = [
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]
    fai    = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    BCFTOOLS_REHEADER ( input, fai )
}
