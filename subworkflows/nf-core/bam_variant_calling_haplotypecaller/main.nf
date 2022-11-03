include { GATK4_HAPLOTYPECALLER  } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { BEDTOOLS_SPLIT         } from '../../../modules/nf-core/bedtools/split/main'
include { BCFTOOLS_MERGE         } from '../../../modules/nf-core/bcftools/merge/main'
include { TABIX_TABIX            } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_HAPLOTYPECALLER {

    take:
    ch_bam              // channel: [mandatory] [meta, bam, bai, bed]
    ch_scatter          // channel: [optional]  [meta, scatter_count]
    ch_fasta            // channel: [mandatory] [fasta]
    ch_fasta_fai        // channel: [mandatory] [fasta_fai]
    ch_dict             // channel: [mandatory] [dict]
    ch_dragstr_models   // channel: [optional]  [meta, dragstr_model]
    ch_dbsnp            // channel: [optional]  [dbsnp]
    ch_dbsnp_tbi        // channel: [optional]  [dbsnp_tbi]

    main:

    ch_versions = Channel.empty()

    //
    // Optional: Scatter the BED files
    //

    ch_input = ch_bam.combine(ch_scatter, by:0)
                    .branch({ meta, bam, bai, bed, scatter_count=[] ->
                        new_meta = meta.clone()
                        new_meta.subwf_scatter_count = scatter_count ?: 1
                        multiple: scatter_count > 1
                            return [ new_meta, bam, bai, bed ]
                        single: scatter_count <= 1 || scatter_count == []
                            return [ new_meta, bam, bai, bed ]
                    })

    ch_input.multiple.view({"multiple: $it"})
    ch_input.single.view({"single: $it"})

    ch_beds =

    BEDTOOLS_SPLIT (
        ch_input.multiple.map({ meta, bam, bai, bed -> [ meta, bed ]}),
        scatter_count
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions.first())

    ch_scattered_bams = ch_bam.combine(BEDTOOLS_SPLIT.out.beds.transpose(), by:0)
                            .map({ meta, bam, bai, full_bed, scattered_bed ->
                                new_meta = meta.clone()
                                new_meta.id = scattered_bed.baseName
                                new_meta.sample = meta.id
                                [ new_meta, bam, bai, scattered_bed ]
                            })

    //
    // Call variants using HaplotypeCaller
    //

    if (ch_dragstr_models){
        ch_bam_dragstr = ch_scattered_bams.combine(ch_dragstr_models, by: 0)
    } else {
        ch_bam_dragstr = ch_scattered_bams
    }

    ch_haplotypecaller_input = ch_bam_dragstr.map({ meta, bam, bai, bed, dragstr_model=[] ->
                                                    [ meta, bam, bai, bed, dragstr_model ]
                                                })

    GATK4_HAPLOTYPECALLER (
        ch_haplotypecaller_input,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_dbsnp,
        ch_dbsnp_tbi
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

    //
    // Optional: Merge scattered VCFs
    //

    ch_merge_input = GATK4_HAPLOTYPECALLER.out.vcf
                        .combine(GATK4_HAPLOTYPECALLER.out.tbi, by:0)
                        .map({ meta, vcf, tbi ->
                            new_meta = [:]
                            new_meta.id = meta.sample
                            [ groupKey(new_meta), vcf, tbi ]
                        })
                        .groupTuple()

    BCFTOOLS_MERGE (
        ch_merge_input,
        [],
        ch_fasta,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions.first())

    TABIX_TABIX (
        BCFTOOLS_MERGE.out.merged_variants
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    ch_called_vcfs      = BCFTOOLS_MERGE.out.merged_variants
    ch_called_vcfs_tbi  = TABIX_TABIX.out.tbi

    ch_called_vcfs      = GATK4_HAPLOTYPECALLER.out.vcf
    ch_called_vcfs_tbi  = GATK4_HAPLOTYPECALLER.out.tbi

    emit:
    vcf      = ch_called_vcfs              // channel: [ meta, vcf ]
    tbi      = ch_called_vcfs_tbi          // channel: [ meta, tbi ]

    versions = ch_versions                 // channel: [ versions.yml ]
}

