include { GATK4_HAPLOTYPECALLER  } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GAWK                   } from '../../../modules/nf-core/gawk/main'
include { BEDTOOLS_SPLIT         } from '../../../modules/nf-core/bedtools/split/main'
include { BCFTOOLS_MERGE         } from '../../../modules/nf-core/bcftools/merge/main'
include { TABIX_TABIX            } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_HAPLOTYPECALLER {

    // IMPORTANT: Add the configuration in nextflow.config file to your modules.config

    take:
    ch_bam              // channel: [mandatory] [meta, bam, bai, bed]
    ch_scatter          // channel: [optional]  [meta, scatter_count]
    fasta               // value:   [mandatory] fasta
    fasta_fai           // value:   [mandatory] fasta_fai
    dict                // value:   [mandatory] dict
    ch_dragstr_models   // channel: [optional]  [meta, dragstr_model]
    dbsnp               // value:   [optional]  dbsnp
    dbsnp_tbi           // value:   [optional]  dbsnp_tbi

    main:

    ch_versions = Channel.empty()

    //
    // Optional: Scatter the BED files
    //

    ch_bam
        .join(ch_dragstr_models ?: Channel.empty(), remainder: true)
        .join(ch_scatter ?: Channel.empty(), remainder: true)
        .branch(
            { meta, bam, bai, bed, dragstr_model, scatter_count ->
                new_meta = meta.clone()
                new_meta.subwf_scatter_count = scatter_count ?: 1

                dragstr = dragstr_model ?: []

                multiple_no_bed: !bed && new_meta.subwf_scatter_count > 1
                    return [ new_meta, bam, bai, dragstr ]
                multiple: new_meta.subwf_scatter_count > 1
                    return [ new_meta, bam, bai, bed, dragstr ]
                single: new_meta.subwf_scatter_count <= 1
                    return [ new_meta, bam, bai, bed, dragstr ]
            }
        )
        .set { ch_input }

    GAWK (
        [[id:'fasta_fai'], fasta_fai],
        []
    )
    ch_versions = ch_versions.mix(GAWK.out.versions)

    ch_input.multiple_no_bed
        .combine(
            GAWK.out.output
                .map(
                    { meta, bed ->
                        [ bed ]
                    }
                )
        )
        .map(
            { meta, bam, bai, dragstr, bed ->
                [ meta, bam, bai, bed, dragstr ]
            }
        )
        .mix(ch_input.multiple)
        .set { ch_to_split }

    BEDTOOLS_SPLIT (
        ch_to_split.map({ meta, bam, bai, bed, dragstr_model -> [ meta, bed ]}),
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions.first())

    ch_to_split.combine(BEDTOOLS_SPLIT.out.beds.transpose(), by:0)
        .map(
            { meta, bam, bai, full_bed, dragstr_model, scattered_bed ->
                new_meta = meta.clone()

                // This logic checks if a scatter has been done on the BED files (to prevent errors with BCFTOOLS_MERGE)
                new_meta.id = full_bed.text == scattered_bed.text ? meta.id : scattered_bed.baseName
                new_meta.subwf_scatter_count = full_bed.text == scattered_bed.text ? 1 : meta.subwf_scatter_count
                new_meta.subwf_sample = meta.id

                [ new_meta, bam, bai, scattered_bed, dragstr_model ]
            }
        )
        .mix(ch_input.single)
        .set { ch_scattered_bams }

    //
    // Call variants using HaplotypeCaller
    //

    GATK4_HAPLOTYPECALLER (
        ch_scattered_bams,
        fasta,
        fasta_fai,
        dict,
        dbsnp,
        dbsnp_tbi
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

    //
    // Optional: Merge scattered VCFs
    //

    GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi)
        .branch(
            { meta, vcf, tbi ->
                new_meta = meta.clone()
                new_meta.remove('subwf_scatter_count')
                new_meta.remove('subwf_sample')
                new_meta.id = meta.subwf_sample

                multiple: meta.subwf_scatter_count > 1
                    [ groupKey(new_meta, meta.subwf_scatter_count), vcf, tbi ]
                single: meta.subwf_scatter_count <= 1
                    [ groupKey(new_meta, meta.subwf_scatter_count), vcf, tbi ]
            }
        )
        .set { ch_caller_output }

    BCFTOOLS_MERGE (
        ch_caller_output.multiple.groupTuple(),
        [],
        fasta,
        fasta_fai
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions.first())

    TABIX_TABIX (
        BCFTOOLS_MERGE.out.merged_variants
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    BCFTOOLS_MERGE.out.merged_variants
        .mix(ch_caller_output.single.map({meta, vcf, tbi -> [meta, vcf]}))
        .set { ch_called_vcfs }

    TABIX_TABIX.out.tbi
        .mix(ch_caller_output.single.map({meta, vcf, tbi -> [meta, tbi]}))
        .set { ch_called_vcfs_tbi }

    emit:
    vcf      = ch_called_vcfs              // channel: [ meta, vcf ]
    tbi      = ch_called_vcfs_tbi          // channel: [ meta, tbi ]

    versions = ch_versions                 // channel: [ versions.yml ]
}

