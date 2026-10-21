// Run GATK4_MUTECT2 in tumor normal mode, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls

include { GATK4_CALCULATECONTAMINATION                                } from '../../../modules/nf-core/gatk4/calculatecontamination'
include { GATK4_FILTERMUTECTCALLS                                     } from '../../../modules/nf-core/gatk4/filtermutectcalls'
include { GATK4_GETPILEUPSUMMARIES as GATK4_GETPILEUPSUMMARIES_NORMAL } from '../../../modules/nf-core/gatk4/getpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES as GATK4_GETPILEUPSUMMARIES_TUMOR  } from '../../../modules/nf-core/gatk4/getpileupsummaries'
include { GATK4_LEARNREADORIENTATIONMODEL                             } from '../../../modules/nf-core/gatk4/learnreadorientationmodel'
include { GATK4_MUTECT2                                               } from '../../../modules/nf-core/gatk4/mutect2'

workflow BAM_SOMATIC_VARIANT_CALLING_GATK {
    take:
    ch_input_tumor // channel: [ val(meta), path(input), path(input_index)]
    ch_input_normal // channel: [ val(meta), path(input), path(input_index)]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai), path(gzi) ]
    ch_dict // channel: [ val(meta), path(dict) ]
    ch_alleles // channel: /path/to/alleles
    ch_alleles_tbi // channel: /path/to/alleles/index
    ch_germline_resource // channel: /path/to/germline/resource
    ch_germline_resource_tbi // channel: /path/to/germline/index
    ch_panel_of_normals // channel: /path/to/panel/of/normals
    ch_panel_of_normals_tbi // channel: /path/to/panel/of/normals/index
    ch_interval_file // channel: /path/to/interval/file
    val_tumor_normal_analysis // boolean: if true, perform tumor-normal analysis

    main:

    if (val_tumor_normal_analysis) {

        ch_mutect2_in = ch_input_tumor
            .combine(ch_input_normal)
            .map { meta_tumor, bam_tumor, bai_tumor, meta_normal, bam_normal, bai_normal -> 
            [[id: "${meta_tumor.id}_${meta_normal.id}" ], [bam_tumor, bam_normal], [bai_tumor, bai_normal]]
            }
            .combine(ch_interval_file)

        // Ensure that the panel of normals is not used for tumor-normal analysis, even if provided
        ch_panel_of_normals = []
        ch_panel_of_normals_tbi = []
    }
    else {
        ch_mutect2_in = ch_input_tumor
            .combine(ch_interval_file)
    }

    GATK4_MUTECT2(
        ch_mutect2_in,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_alleles,
        ch_alleles_tbi,
        ch_germline_resource,
        ch_germline_resource_tbi,
        ch_panel_of_normals,
        ch_panel_of_normals_tbi,
    )

    // Generate artifactpriors using learnreadorientationmodel on the f1r2 output of GATK4_MUTECT2
    GATK4_LEARNREADORIENTATIONMODEL(GATK4_MUTECT2.out.f1r2)

    // Generate pileup summary tables using getepileupsummaries
    ch_getpileupsummaries_tumor_in = ch_input_tumor
            .combine(ch_interval_file)

    GATK4_GETPILEUPSUMMARIES_TUMOR(
        ch_getpileupsummaries_tumor_in,
        ch_fasta,
        ch_fai.map { meta, fai, _gzi -> [meta, fai] },
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_tbi,
    )

    ch_pileup_table_tumor = GATK4_GETPILEUPSUMMARIES_TUMOR.out.table

    if (val_tumor_normal_analysis) {

        ch_getpileupsummaries_normal_in = ch_input_normal
            .combine(ch_interval_file)

        GATK4_GETPILEUPSUMMARIES_NORMAL(
            ch_getpileupsummaries_normal_in,
            ch_fasta,
            ch_fai.map { meta, fai, _gzi -> [meta, fai] },
            ch_dict,
            ch_germline_resource,
            ch_germline_resource_tbi,
        )

        ch_pileup_table_normal = GATK4_GETPILEUPSUMMARIES_NORMAL.out.table

        // Ensure the correct order of the files in the input channel for GATK4_CALCULATECONTAMINATION
        ch_calccon_in = ch_pileup_table_tumor
            .combine(ch_pileup_table_normal)
            .map { meta_tumor, table_tumor, meta_normal, table_normal -> 
                [[id: "${meta_tumor.id}_${meta_normal.id}"], table_tumor, table_normal]
            }
    }
    else {
        ch_pileup_table_normal = channel.empty()

        ch_calccon_in = ch_pileup_table_tumor
            .map { meta_tumor, table_tumor -> 
                [meta_tumor, table_tumor, []]
            }
    }

    // Contamination and segmentation tables created using calculatecontamination on the pileup summary table
    GATK4_CALCULATECONTAMINATION(ch_calccon_in)

    // GATK4_MUTECT2 calls filtered by filtermutectcalls using the artifactpriors, contamination and segmentation tables
    ch_filtermutect_in = GATK4_MUTECT2.out.vcf
        .join(GATK4_MUTECT2.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(GATK4_MUTECT2.out.stats, failOnDuplicate: true, failOnMismatch: true)
        .join(GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior, failOnDuplicate: true, failOnMismatch: true)
        .join(GATK4_CALCULATECONTAMINATION.out.segmentation, failOnDuplicate: true, failOnMismatch: true)
        .join(GATK4_CALCULATECONTAMINATION.out.contamination, failOnDuplicate: true, failOnMismatch: true) 
        .map{ meta, vcf, tbi, stats, artifactprior, segmentation, contamination ->
            [   
                meta,
                vcf,
                tbi,
                stats,
                artifactprior,
                segmentation,
                contamination,
                [], // Since the contamination table is provided, the value estimate is not relevant
            ]
        }

    GATK4_FILTERMUTECTCALLS(
        ch_filtermutect_in,
        ch_fasta,
        ch_fai.map { meta, fai, _gzi -> [meta, fai] },
        ch_dict,
    )

    emit:
    artifact_priors     = GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior // channel: [ val(meta), path(artifactprior) ]
    contamination_table = GATK4_CALCULATECONTAMINATION.out.contamination // channel: [ val(meta), path(table) ]
    filtered_stats      = GATK4_FILTERMUTECTCALLS.out.stats // channel: [ val(meta), path(stats) ]
    filtered_tbi        = GATK4_FILTERMUTECTCALLS.out.tbi // channel: [ val(meta), path(tbi) ]
    filtered_vcf        = GATK4_FILTERMUTECTCALLS.out.vcf // channel: [ val(meta), path(vcf) ]
    mutect2_f1r2        = GATK4_MUTECT2.out.f1r2 // channel: [ val(meta), path(f1r2) ]
    mutect2_stats       = GATK4_MUTECT2.out.stats // channel: [ val(meta), path(stats) ]
    mutect2_tbi         = GATK4_MUTECT2.out.tbi // channel: [ val(meta), path(tbi) ]
    mutect2_vcf         = GATK4_MUTECT2.out.vcf // channel: [ val(meta), path(vcf) ]
    pileup_table_normal = ch_pileup_table_normal // channel: [ val(meta), path(table) ]
    pileup_table_tumor  = ch_pileup_table_tumor // channel: [ val(meta), path(table) ]
    segmentation_table  = GATK4_CALCULATECONTAMINATION.out.segmentation // channel: [ val(meta), path(table) ]
}