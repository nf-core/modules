#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOMALIER_EXTRACT } from '../../../../../modules/nf-core/somalier/extract/main.nf'

workflow test_somalier_extract {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam_bai'], checkIfExists: true)
    ]

    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fasta_fai   = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)

    sites       = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/sites_chr21.hg38.vcf.gz", checkIfExists: true)

    SOMALIER_EXTRACT ( input, fasta, fasta_fai, sites )
}
