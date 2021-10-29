#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MULTIVCFANALYZER } from '../../../modules/multivcfanalyzer/main.nf' addParams( options: [:] )

workflow test_multivcfanalyzer {
    
    vcfs = [ file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true),
             file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true) ]

    reference = 
    allele_freq = T
    minq = 30
    min_cov = 5
    min_freq_homo = 0.9
    min_freq_hetero = 0.1

    MULTIVCFANALYZER ( vcfs, reference, allele_freq, minq, min_cov, min_freq_homo, min_freq_hetero )
}
