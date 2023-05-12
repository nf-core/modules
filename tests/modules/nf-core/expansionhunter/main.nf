#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EXPANSIONHUNTER } from '../../../../modules/nf-core/expansionhunter/main.nf'

workflow test_expansionhunter {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
            ]
    fasta = [[id:'fasta'],file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[id:'fasta_fai'],file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]
    variant_catalog = [[id:'catalogue'],file(params.test_data['homo_sapiens']['genome']['repeat_expansions'], checkIfExists: true)]

    EXPANSIONHUNTER (
      input,
      fasta,
      fasta_fai,
      variant_catalog
    )
}
