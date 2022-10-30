#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_UNIFIEDGENOTYPER } from '../../../../modules/nf-core/gatk/unifiedgenotyper/main.nf'
include { GUNZIP                } from '../../../../modules/nf-core/gunzip/main.nf'
include { MULTIVCFANALYZER      } from '../../../../modules/nf-core/multivcfanalyzer/main.nf'

workflow test_multivcfanalyzer {

    input     = Channel.of([ [ id:'test' ], // meta map
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                ],
                [ [ id:'test2' ], // meta map
                  file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
                  file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true),
                ],
                )
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK_UNIFIEDGENOTYPER ( input, fasta, fai, dict, [], [], [], [])

    mva_vcf = GUNZIP ( GATK_UNIFIEDGENOTYPER.out.vcf ).gunzip
        .map{it[1]}
        .collect()

    snpeff_results          = []
    gff                     = []
    allele_freqs            = true
    genotype_quality        = 30
    coverage                = 5
    homozygous_freq         = 0.8
    heterozygous_freq       = 0.2
    gff_exclude             = []

    MULTIVCFANALYZER ( mva_vcf, fasta, snpeff_results, gff, allele_freqs, genotype_quality, coverage, homozygous_freq, heterozygous_freq, gff_exclude )
}
