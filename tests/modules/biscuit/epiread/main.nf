#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISCUIT_INDEX   } from '../../../../modules/biscuit/index/main.nf'
include { BISCUIT_EPIREAD } from '../../../../modules/biscuit/epiread/main.nf'

workflow test_biscuit_epiread_nosnp {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam_bai'], checkIfExists: true),
        [] //SNP BED file
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BISCUIT_INDEX( fasta )
    BISCUIT_EPIREAD ( input, BISCUIT_INDEX.out.index )
}

workflow test_biscuit_epiread_snp {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam_bai'], checkIfExists: true),
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/biscuit/test-snp.bed')
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BISCUIT_INDEX( fasta )
    BISCUIT_EPIREAD ( input, BISCUIT_INDEX.out.index )
}

workflow test_biscuit_epiread_snp_decompress {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam_bai'], checkIfExists: true),
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/biscuit/test-snp.bed.gz')
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BISCUIT_INDEX( fasta )
    BISCUIT_EPIREAD ( input, BISCUIT_INDEX.out.index )
}
