#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ATLAS_PMD } from '../../../../../modules/nf-core/atlas/pmd/main.nf'

workflow test_atlas_pmd {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        []
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    ATLAS_PMD ( input, fasta, fai )
}
