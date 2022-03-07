#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ASCAT as ASCAT_SIMPLE} from '../../../modules/ascat/main.nf'
include { ASCAT as ASCAT_PLOIDY_AND_PURITY} from '../../../modules/ascat/main.nf'
include { ASCAT as ASCAT_CRAM} from '../../../modules/ascat/main.nf'




workflow test_ascat {


    input = [
        [ id:'test', single_end:false ], // meta map
        file("/home/ec2-user/input_files/bams/HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam", checkIfExists: true),
        file("/home/ec2-user/input_files/bams/HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam.bai", checkIfExists: true),
        file("/home/ec2-user/input_files/bams/HG00155.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam", checkIfExists: true),
        file("/home/ec2-user/input_files/bams/HG00155.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam.bai", checkIfExists: true)
    ]

    ASCAT_SIMPLE ( input , "/home/ec2-user/input_files/allele_files", "/home/ec2-user/input_files/loci_files")
}








workflow test_ascat_with_ploidy_and_purity {  
   input = [
        [ id:'test', single_end:false ], // meta map
        file("/home/ec2-user/input_files/bams/HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam", checkIfExists: true),
        file("/home/ec2-user/input_files/bams/HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam.bai", checkIfExists: true),
        file("/home/ec2-user/input_files/bams/test2.bam", checkIfExists: true),
        file("/home/ec2-user/input_files/bams/test2.bam.bai", checkIfExists: true)
    ]

    ASCAT_PLOIDY_AND_PURITY ( input , "/home/ec2-user/input_files/allele_files", "/home/ec2-user/input_files/loci_files")
}



workflow test_ascat_with_crams {


    input = [
        [ id:'test', single_end:false ], // meta map
        file("/home/ec2-user/input_files/crams/HG00145.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram", checkIfExists: true),
        file("/home/ec2-user/input_files/crams/HG00145.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram.crai", checkIfExists: true),
        file("/home/ec2-user/input_files/crams/duplicate_test.cram", checkIfExists: true),
        file("/home/ec2-user/input_files/crams/duplicate_test.cram.crai", checkIfExists: true)
    ]

    ASCAT_CRAM ( input , "/home/ec2-user/input_files/allele_files", "/home/ec2-user/input_files/loci_files")
}



