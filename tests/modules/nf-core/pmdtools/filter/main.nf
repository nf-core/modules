#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PMDTOOLS_FILTER } from '../../../../modules/pmdtools/filter/main.nf'

workflow test_pmdtools_filter {
    
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]]
    threshold = 3
    reference = [ file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]
    PMDTOOLS_FILTER ( input, threshold, reference )
}
