#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_COLLECTREADCOUNTS                                                               } from '../../../../../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY as GATK4_DETERMINEGERMLINECONTIGPLOIDY_COHORT     } from '../../../../../modules/nf-core/gatk4/determinegermlinecontigploidy/main.nf'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY as GATK4_DETERMINEGERMLINECONTIGPLOIDY_CASE       } from '../../../../../modules/nf-core/gatk4/determinegermlinecontigploidy/main.nf'
include { GATK4_GERMLINECNVCALLER as GATK4_GERMLINECNVCALLER_COHORT                             } from '../../../../../modules/nf-core/gatk4/germlinecnvcaller/main.nf'
include { GATK4_GERMLINECNVCALLER as GATK4_GERMLINECNVCALLER_CASE                               } from '../../../../../modules/nf-core/gatk4/germlinecnvcaller/main.nf'
include { GATK4_POSTPROCESSGERMLINECNVCALLS                                                     } from '../../../../../modules/nf-core/gatk4/postprocessgermlinecnvcalls/main.nf'

workflow test_gatk4_postprocessgermlinecnvcalls {
    bed = Channel.of(file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true))
    priors =  file(params.test_data['homo_sapiens']['illumina']['contig_ploidy_priors_table'], checkIfExists: true)
    inputs = Channel.of([
        [ id:'test1' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ],
    [
        [ id:'test2' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ])
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    ploidy_calls = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_ploidy_calls'], checkIfExists: true) ])
    cnv_calls = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_germline_cnv_calls'], checkIfExists: true) ])
    cnv_model = Channel.of(file(params.test_data['homo_sapiens']['genome']['genome_germline_cnv_model'], checkIfExists: true))

    GATK4_POSTPROCESSGERMLINECNVCALLS ( ploidy_calls, cnv_model, cnv_calls.collect{ it[1] } )
    GATK4_POSTPROCESSGERMLINECNVCALLS.out.intervals
    GATK4_POSTPROCESSGERMLINECNVCALLS.out.segments
    GATK4_POSTPROCESSGERMLINECNVCALLS.out.denoised
}
