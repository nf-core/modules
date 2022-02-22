#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ASCAT } from '../../../modules/ascat/main.nf'

workflow test_ascat {
    
//    input = [ 
//        [ id:'test', single_end:false ], // meta map
//        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
//        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
//        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
//        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)
//    ]


    input = [
        [ id:'test', single_end:false ], // meta map
        file("/home/ec2-user/input_files/bams/HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam", checkIfExists: true),
        file("/home/ec2-user/input_files/bams/HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam.bai", checkIfExists: true),
        file("/home/ec2-user/input_files/bams/test2.bam", checkIfExists: true),
        file("/home/ec2-user/input_files/bams/test2.bam.bai", checkIfExists: true)
    ]

    ASCAT ( input , "/home/ec2-user/input_files/allele_files", "/home/ec2-user/input_files/loci_files")
}
