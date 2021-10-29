#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ATAQV_ATAQV } from '../../../../modules/ataqv/ataqv/main.nf' addParams( options: [:] )
include { ATAQV_ATAQV as ATAQV_ATAQV_PROBLEM_READS} from '../../../../modules/ataqv/ataqv/main.nf' addParams( options: ['args': '--log-problematic-reads'] )

workflow test_ataqv_ataqv {
    
    meta_map = [ id:'test', single_end:false ] 
    bam_file = file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    bai_file = []
    peak_file = []
    organism = 'human'
    tss_file = []
    excl_regs_file = []
    autosom_ref_file = []

    input_samp_spec = [ meta_map,
                        bam_file,
                        bai_file, 
                        peak_file ]

    ATAQV_ATAQV ( input_samp_spec, organism, tss_file, excl_regs_file, autosom_ref_file )
}

workflow test_ataqv_ataqv_problem_reads {
    
    meta_map = [ id:'test', single_end:false ] 
    bam_file = file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    bai_file = []
    peak_file = []
    organism = 'human'
    tss_file = []
    excl_regs_file = []
    autosom_ref_file = []

    input_samp_spec = [ meta_map,
                        bam_file,
                        bai_file, 
                        peak_file ]

    ATAQV_ATAQV_PROBLEM_READS ( input_samp_spec, organism, tss_file, excl_regs_file, autosom_ref_file )
}

workflow test_ataqv_ataqv_peak {
    
    meta_map = [ id:'test', single_end:false ] 
    bam_file = file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    bai_file = []
    peak_file = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    organism = 'human'
    tss_file = []
    excl_regs_file = []
    autosom_ref_file = []

    input_samp_spec = [ meta_map,
                        bam_file,
                        bai_file, 
                        peak_file ]

    ATAQV_ATAQV ( input_samp_spec, organism, tss_file, excl_regs_file, autosom_ref_file )
}

workflow test_ataqv_ataqv_tss {
    
    meta_map = [ id:'test', single_end:false ] 
    bam_file = file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    bai_file = file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    peak_file = []
    organism = 'human'
    tss_file = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    excl_regs_file = []
    autosom_ref_file = []

    input_samp_spec = [ meta_map,
                        bam_file,
                        bai_file, 
                        peak_file ]

    ATAQV_ATAQV ( input_samp_spec, organism, tss_file, excl_regs_file, autosom_ref_file )
}

workflow test_ataqv_ataqv_excluded_regs {
    
    meta_map = [ id:'test', single_end:false ] 
    bam_file = file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    bai_file = file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    peak_file = []
    organism = 'human'
    tss_file = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    excl_regs_file = file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
    autosom_ref_file = []

    input_samp_spec = [ meta_map,
                        bam_file,
                        bai_file, 
                        peak_file ]

    ATAQV_ATAQV ( input_samp_spec, organism, tss_file, excl_regs_file, autosom_ref_file )
}
