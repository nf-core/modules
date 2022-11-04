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

    ch_bam
        .join(ch_dragstr_models ?: Channel.empty(), remainder: true)
        .join(ch_scatter ?: Channel.empty(), remainder: true)
        .branch( 
            { meta, bam, bai, bed, dragstr_model=[], scatter_count=[] ->
                new_meta = meta.clone()
                new_meta.subwf_scatter_count = scatter_count ?: 1

                dragstr = dragstr_model ?: []

                multiple: new_meta.subwf_scatter_count > 1
                    return [ new_meta, bam, bai, bed, dragstr ]
                single: new_meta.subwf_scatter_count <= 1
                    return [ new_meta, bam, bai, bed, dragstr ]
            }
        )
        .set { ch_input }

    BEDTOOLS_SPLIT (
        ch_input.multiple.map({ meta, bam, bai, bed, dragstr_model -> [ meta, bed ]}),
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions.first())

    ch_scattered_bams = ch_input.multiple.combine(BEDTOOLS_SPLIT.out.beds.transpose(), by:0)
                            .map({ meta, bam, bai, full_bed, dragstr_model, scattered_bed ->
                                new_meta = meta.clone()
                                new_meta.id = scattered_bed.baseName
                                new_meta.subwf_sample = meta.id
                                [ new_meta, bam, bai, scattered_bed, dragstr_model ]
                            })
                            .mix(ch_input.single)

    //
    // Call variants using HaplotypeCaller
    //

    GATK4_HAPLOTYPECALLER (
        ch_scattered_bams,
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

    ch_caller_output = GATK4_HAPLOTYPECALLER.out.vcf
                        .join(GATK4_HAPLOTYPECALLER.out.tbi)
                        .branch({ meta, vcf, tbi ->
                            new_meta = meta.clone()
                            new_meta.remove('subwf_scatter_count')
                            new_meta.remove('subwf_sample')
                            new_meta.id = meta.subwf_sample

                            multiple: meta.subwf_scatter_count > 1
                                [ groupKey(new_meta, meta.subwf_scatter_count), vcf, tbi ]
                            single: meta.subwf_scatter_count <= 1
                                [ groupKey(new_meta, meta.subwf_scatter_count), vcf, tbi ]
                        })

    BCFTOOLS_MERGE (
        ch_caller_output.multiple.groupTuple(),
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
                            .mix(ch_caller_output.single.map({meta, vcf, tbi -> [meta, vcf]}))
    ch_called_vcfs_tbi  = TABIX_TABIX.out.tbi
                            .mix(ch_caller_output.single.map({meta, vcf, tbi -> [meta, tbi]}))

    emit:
    vcf      = ch_called_vcfs              // channel: [ meta, vcf ]
    tbi      = ch_called_vcfs_tbi          // channel: [ meta, tbi ]

    versions = ch_versions                 // channel: [ versions.yml ]
}

