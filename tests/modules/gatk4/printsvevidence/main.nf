#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_PRINTSVEVIDENCE as GATK4_PRINTSVEVIDENCE_PE } from '../../../../modules/gatk4/printsvevidence/main.nf'
include { GATK4_PRINTSVEVIDENCE as GATK4_PRINTSVEVIDENCE_SR } from '../../../../modules/gatk4/printsvevidence/main.nf'
include { GATK4_COLLECTSVEVIDENCE                           } from '../../../../modules/gatk4/collectsvevidence/main.nf'

workflow test_gatk4_printsvevidence {
    
    input = Channel.of(
        [
            [ id:'normal', single_end:false ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
            [],
            []
        ],
        [
            [ id:'tumour', single_end:false ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
            [],
            []
        ]
    )

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)


    GATK4_COLLECTSVEVIDENCE(
        input, 
        fasta, 
        fasta_fai, 
        dict
    )

    pe_files = GATK4_COLLECTSVEVIDENCE.out.paired_end_evidence
                    .combine(GATK4_COLLECTSVEVIDENCE.out.paired_end_evidence_index, by:0)
    sr_files = GATK4_COLLECTSVEVIDENCE.out.split_read_evidence
                    .combine(GATK4_COLLECTSVEVIDENCE.out.split_read_evidence_index, by:0)

    GATK4_PRINTSVEVIDENCE_PE( 
        pe_files,
        [],
        fasta,
        fasta_fai,
        dict
    )

    GATK4_PRINTSVEVIDENCE_SR( 
        sr_files,
        [],
        fasta,
        fasta_fai,
        dict
    )

}

workflow test_gatk4_printsvevidence_bed_no_fasta {
    
    input = Channel.of(
        [
            [ id:'normal', single_end:false ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
            [],
            []
        ],
        [
            [ id:'tumour', single_end:false ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
            [],
            []
        ]
    )

    fasta = []
    fasta_fai = []
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    bed = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)

    GATK4_COLLECTSVEVIDENCE(
        input, 
        fasta, 
        fasta_fai, 
        dict
    )

    pe_files = GATK4_COLLECTSVEVIDENCE.out.paired_end_evidence
                    .combine(GATK4_COLLECTSVEVIDENCE.out.paired_end_evidence_index, by:0)
    sr_files = GATK4_COLLECTSVEVIDENCE.out.split_read_evidence
                    .combine(GATK4_COLLECTSVEVIDENCE.out.split_read_evidence_index, by:0)

    GATK4_PRINTSVEVIDENCE_PE( 
        pe_files,
        bed,
        fasta,
        fasta_fai,
        dict
    )

    GATK4_PRINTSVEVIDENCE_SR( 
        sr_files,
        bed,
        fasta,
        fasta_fai,
        dict
    )

}
