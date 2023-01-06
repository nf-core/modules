#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
params.enable_conda = true

include { GRAPHTYPER_GENOTYPE } from '../../../../../modules/nf-core/graphtyper/genotype/main.nf'

workflow test_graphtyper_genotype {
    
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]
    ]
    reference = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ref_index = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    region = [ 'MT192765.1' ]

    GRAPHTYPER_GENOTYPE ( input, reference, ref_index, region )
}
