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

    main:
    ch_versions    = Channel.empty()

    // Merge bed files to create a single interval list
    // BEDTOOLS_MULTIINTER( [meta, bed], chrom_sizes )
    BEDTOOLS_MULTIINTER( ch_beds, [] )

    // Split merged bed into individual intervals for genotyping
    // BED_SCATTER_BEDTOOLS( [meta, bed, scatter_count] )
    BED_SCATTER()

    // Generate (multi)sample GenomicsDB
    // GATK4_GENOMICSDBIMPORT([ meta, vcf, vcf_index, interval_file, interval_value, workspace ], run_intlist, run_updatewspace, input_map)
    GATK4_GENOMICSDBIMPORT(
        ch_gvcfs.map{meta, gvcf, tbi -> [meta, gvcf, tbi, [], [], []]},
        false,
        false,
        []
    )

    // GenotypeGVCFs
    GATK4_GENOTYPEGVCFS (genotype_input, fasta, fai, dict, dbsnp, dbsnp_tbi)

    // Merge split GenotypeGVCFs VCFs
    GATK4_MERGEVCFS(
        GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{ meta, gdb -> [meta, gdb, [], [], []]},
        ch_dict
    )

    // Generate VCF index
    TABIX()

    //
    //Map input for GenomicsDBImport.
    //rename based on num_intervals, group all samples by their interval_name/interval_file and restructure for channel
    //group by 0,3 to avoid a list of metas and make sure that any intervals
    //

    // Generate genomicsdb for each interval in the sample
    gendb_input = input.map{
        meta, gvcf, tbi, intervals->
            new_meta = [
                        id:             "joint_variant_calling",
                        intervals_name: meta.intervals_name,
                        num_intervals:  meta.num_intervals
                    ]

            [ new_meta, gvcf, tbi, intervals ]

        }.groupTuple(by:[0,3]).map{ new_meta, gvcf, tbi, intervals ->
            [ new_meta, gvcf, tbi, intervals, [], [] ]
        }

    //
    //Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport.
    //
    GATK4_GENOMICSDBIMPORT ( gendb_input, false, false, false )

    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{
        meta, genomicsdb ->
            [meta, genomicsdb, [], [], []]
        }

    //
    //Joint genotyping performed using GenotypeGVCFs
    //Sort vcfs called by interval within each VCF
    //

    GATK4_GENOTYPEGVCFS (genotype_input, fasta, fai, dict, dbsnp, dbsnp_tbi)

    BCFTOOLS_SORT(GATK4_GENOTYPEGVCFS.out.vcf)
    vcfs_sorted_input = BCFTOOLS_SORT.out.vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    vcfs_sorted_input_no_intervals =  vcfs_sorted_input.no_intervals.map{ meta , vcf ->

                            [[  id:             "joint_variant_calling",
                                num_intervals:  meta.num_intervals,
                                patient:        "all_samples",
                                variantcaller:  "haplotypecaller"
                            ] , vcf ]
    }

    // Index vcf files if no scatter/gather by intervals
    TABIX(vcfs_sorted_input_no_intervals)

    //Merge scatter/gather vcfs & index
    //Rework meta for variantscalled.csv and annotation tools
    MERGE_GENOTYPEGVCFS(vcfs_sorted_input.intervals.map{meta, vcf ->
                                [
                                        [
                                            id: "joint_variant_calling",
                                            num_intervals: meta.num_intervals,
                                            patient: "all_samples",
                                            variantcaller: "haplotypecaller",
                                        ],
                                vcf]
                            }.groupTuple()
                            ,dict)

    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)




    emit:
    versions       = ch_versions                                                                                            // channel: [ versions.yml ]
    genotype_vcf   = Channel.empty().mix( vcfs_sorted_input_no_intervals, MERGE_GENOTYPEGVCFS.out.vcf)  // channel: [ val(meta), [ vcf ] ]
    genotype_index = Channel.empty().mix( TABIX.out.tbi, MERGE_GENOTYPEGVCFS.out.tbi)                    // channel: [ val(meta), [ tbi ] ]
}

