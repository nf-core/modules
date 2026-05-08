// Run GATK4_MUTECT2 in tumor normal mode, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls

include { GATK4_CALCULATECONTAMINATION                                } from '../../../modules/nf-core/gatk4/calculatecontamination'
include { GATK4_FILTERMUTECTCALLS                                     } from '../../../modules/nf-core/gatk4/filtermutectcalls'
include { GATK4_GETPILEUPSUMMARIES as GATK4_GETPILEUPSUMMARIES_NORMAL } from '../../../modules/nf-core/gatk4/getpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES as GATK4_GETPILEUPSUMMARIES_TUMOR  } from '../../../modules/nf-core/gatk4/getpileupsummaries'
include { GATK4_LEARNREADORIENTATIONMODEL                             } from '../../../modules/nf-core/gatk4/learnreadorientationmodel'
include { GATK4_MUTECT2                                               } from '../../../modules/nf-core/gatk4/mutect2'

workflow BAM_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING_GATK {
    take:
    ch_input // channel: [ val(meta), path(input), path(input_index), val(which_norm) ]
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

    main:
    // Perform variant calling using GATK4_MUTECT2 module in tumor single mode.
    GATK4_MUTECT2(
        ch_input,
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

    // Generate artifactpriors using learnreadorientationmodel on the f1r2 output of GATK4_MUTECT2.
    GATK4_LEARNREADORIENTATIONMODEL(GATK4_MUTECT2.out.f1r2.collect())

    // Generate pileup summary tables using getepileupsummaries
    // Tumor sample should always be passed in as the first input and input list entries of ch_input,
    // to ensure correct file order for calculatecontamination.
    ch_pileup_tumor_input = ch_input
        .combine(ch_interval_file)
        .map { meta, input_file, input_index, _which_norm, intervals ->
            [meta, input_file[0], input_index[0], intervals]
        }

    ch_pileup_normal_input = ch_input
        .combine(ch_interval_file)
        .map { meta, input_file, input_index, _which_norm, intervals ->
            [meta, input_file[1], input_index[1], intervals]
        }

    GATK4_GETPILEUPSUMMARIES_TUMOR(
        ch_pileup_tumor_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_tbi,
    )

    GATK4_GETPILEUPSUMMARIES_NORMAL(
        ch_pileup_normal_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_tbi,
    )

    // Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    ch_pileup_tumor = GATK4_GETPILEUPSUMMARIES_TUMOR.out.table.collect()
    ch_pileup_normal = GATK4_GETPILEUPSUMMARIES_NORMAL.out.table.collect()
    ch_calccon_in = ch_pileup_tumor.join(ch_pileup_normal, failOnDuplicate: true, failOnMismatch: true)

    GATK4_CALCULATECONTAMINATION(ch_calccon_in)

    // GATK4_MUTECT2 calls filtered by filtermutectcalls using the artifactpriors, contamination and segmentation tables.
    ch_vcf = GATK4_MUTECT2.out.vcf.collect()
    ch_tbi = GATK4_MUTECT2.out.tbi.collect()
    ch_stats = GATK4_MUTECT2.out.stats.collect()
    ch_orientation = GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior.collect()
    ch_segment = GATK4_CALCULATECONTAMINATION.out.segmentation.collect()

    // [] is used as a placeholder for optional input to specify the contamination estimate as a value, since the contamination table is used, this is not needed.
    ch_contamination = GATK4_CALCULATECONTAMINATION.out.contamination.map { meta, table -> [meta, table, []] }.collect()

    ch_filtermutect_in = ch_vcf
        .join(ch_tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_stats, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_orientation, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_segment, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_contamination, failOnDuplicate: true, failOnMismatch: true)

    GATK4_FILTERMUTECTCALLS(
        ch_filtermutect_in,
        ch_fasta,
        ch_fai,
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
    pileup_table_normal = GATK4_GETPILEUPSUMMARIES_NORMAL.out.table // channel: [ val(meta), path(table) ]
    pileup_table_tumor  = GATK4_GETPILEUPSUMMARIES_TUMOR.out.table // channel: [ val(meta), path(table) ]
    segmentation_table  = GATK4_CALCULATECONTAMINATION.out.segmentation // channel: [ val(meta), path(table) ]
}
