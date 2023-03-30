#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NGSCHECKMATE_NCM as NGSCHECKMATE_NCM_BAM} from '../../../../../modules/nf-core/ngscheckmate/ncm/main.nf'
include { NGSCHECKMATE_NCM as NGSCHECKMATE_NCM_VCF} from '../../../../../modules/nf-core/ngscheckmate/ncm/main.nf'

include { BEDTOOLS_MAKEWINDOWS } from '../../../../../modules/nf-core/bedtools/makewindows/main.nf'

include { BCFTOOLS_MPILEUP } from '../../../../../modules/nf-core/bcftools/mpileup/main.nf'
include { BCFTOOLS_MPILEUP as BCFTOOLS_MPILEUP_TWO } from '../../../../../modules/nf-core/bcftools/mpileup/main.nf'

workflow test_ngscheckmate_ncm_bam {
    input    = [ file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam'], checkIfExists: true),
                 file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam_bai'], checkIfExists: true),
                 file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                 file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)]

    fasta    = [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    inputBed = [ [ id:'test'],
                 file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)]

    BEDTOOLS_MAKEWINDOWS(inputBed).
    bed.
    map{it[1]}.
    view().
    set{snp_channel}

    NGSCHECKMATE_NCM_BAM(input, snp_channel, fasta)
}

workflow test_ngscheckmate_ncm_vcf {
    fasta    = [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    inputBed = [ [ id:'test'],
                file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)]

    input1   = [ [ id:'test1' ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
                ]

    input2   = [ [ id:'test2' ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
                ]

    BCFTOOLS_MPILEUP ( input1, fasta, false )
    BCFTOOLS_MPILEUP_TWO ( input2, fasta, false )

    BCFTOOLS_MPILEUP_TWO.out.vcf.
        combine( BCFTOOLS_MPILEUP.out.vcf ).
        map { [ it[1], it[3] ] }.
        set { vcf_channel }

    BEDTOOLS_MAKEWINDOWS( inputBed).bed.
        map { it[1] }.
        view().
        set { snp_channel }

    NGSCHECKMATE_NCM_VCF(vcf_channel, snp_channel, fasta)
}


