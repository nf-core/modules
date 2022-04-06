#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DAMAGEPROFILER } from '../../../modules/damageprofiler/main.nf'

workflow test_damageprofiler {

    input        = [ [ id:'test', single_end:false ], // meta map
                   [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam'], checkIfExists: true) ] ]
    fasta        = []
    fai          = []
    species_list = []


    DAMAGEPROFILER ( input, fasta, fai, species_list )
}

workflow test_damageprofiler_reference {

    input        = [ [ id:'test', single_end:false ], // meta map
                   [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam'], checkIfExists: true) ] ]
    fasta        = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai          = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    species_list = []

    DAMAGEPROFILER ( input, fasta, fai, species_list )
}

workflow test_damageprofiler_specieslist {

    input        = [ [ id:'test', single_end:false ], // meta map
                   [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam'], checkIfExists: true) ] ]
    fasta        = []
    fai          = []
    species_list = file(params.test_data['homo_sapiens']['genome']['genome_header'], checkIfExists: true)

    DAMAGEPROFILER ( input, fasta, fai, species_list )
}
