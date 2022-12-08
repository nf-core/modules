#!/usr/bin/env nextflow

include { GATK4_GENOMICSDBIMPORT                                } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS                                   } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATK4_MERGEVCFS as MERGE_GENOTYPEGVCFS                } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { TABIX_TABIX as TABIX                                  } from '../../../modules/nf-core/tabix/tabix/main'
include { BEDTOOLS_MULTIINTER                                   } from '../../../modules/nf-core/bedtools/multiinter/main'

include { BED_SCATTER_BEDTOOLS as BED_SCATTER                   } from '../../subworkflows/nf-core/bed_scatter_bedtools/main'

workflow VCF_JOINT_CALLING_GERMLINE_GATK {

    take:
    ch_gvcfs            // channel: [ val(meta), [ gvcf ], [ gvcf_index ]]
    ch_beds             // channel: [ val(meta), [ bed ]]
    ch_fasta            // channel: [ val(meta), /path/to/reference/fasta]
    ch_fai              // channel: [ val(meta), /path/to/reference/fasta/index]
    ch_dict             // channel: [ val(meta), /path/to/reference/fasta/dict]
    ch_dbsnp            // channel: [ val(meta), /path/to/dbsnp/vcf]
    ch_dbsnp_tbi        // channel: [ val(meta), /path/to/dbsnp/vcf/index]
    scatter_count       // channel: val(scatter_count)

    main:
    ch_versions    = Channel.empty()

    // Merge bed files to create a single interval list
    // BEDTOOLS_MULTIINTER( [meta, bed], chrom_sizes )
    BEDTOOLS_MULTIINTER( ch_beds, [] )
    ch_versions = ch_versions.mix(BEDTOOLS_MULTIINTER.out.versions).first()

    // Split merged bed into individual intervals for genotyping
    // BED_SCATTER_BEDTOOLS( [meta, bed, scatter_count] )
    BED_SCATTER_BEDTOOLS(BEDTOOLS_MULTIINTER.out.bed, scatter_count)
    ch_versions = ch_versions.mix(BED_SCATTER_BEDTOOLS.out.versions).first()

    // Generate (multi)sample GenomicsDB
    // GATK4_GENOMICSDBIMPORT([ meta, vcf, vcf_index, interval_file, interval_value, workspace ], run_intlist, run_updatewspace, input_map)
    GATK4_GENOMICSDBIMPORT(
        ch_gvcfs.map{meta, gvcf, tbi -> [meta, gvcf, tbi, [], [], []]},
        false,
        false,
        []
    )
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions).first()

    // Join gdb with each interval
    ch_gdb_with_interval = GATK4_GENOMICSDBIMPORT.out.genomicsdb.combine(BED_SCATTER_BEDTOOLS.out.scattered_beds)

    // GenotypeGVCFs
    GATK4_GENOTYPEGVCFS (ch_gdb_with_interval, fasta, fai, dict, dbsnp, dbsnp_tbi)
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions).first()

    // Merge split GenotypeGVCFs VCFs
    GATK4_MERGEVCFS(
        GATK4_GENOTYPEGVCFS.out.vcf.map{ meta, gdb -> [meta, gdb, [], [], []]},
        ch_dict
    )
    ch_versions = ch_versions.mix(GATK4_MERGEVCFS.out.versions).first()

    // Sort merged VCF
    BCFTOOLS_SORT(GATK4_MERGEVCFS.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions).first()

    // Generate VCF index
    TABIX(BCFTOOLS_SORT.out.vcf)
    ch_versions = ch_versions.mix(TABIX.out.versions).first()


    emit:
    versions    = ch_versions                                                                                            // channel: [ versions.yml ]
    vcf_tbi     = TABIX.out.vcf_tbi                                                                                      // channel: [ val(meta), vcf, index ]
}

