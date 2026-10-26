//
// Run GATK4_MUTECT2, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls for tumor normal or tumor-only analyses
//

include { GATK4_CALCULATECONTAMINATION    } from '../../../modules/nf-core/gatk4/calculatecontamination'
include { GATK4_FILTERMUTECTCALLS         } from '../../../modules/nf-core/gatk4/filtermutectcalls'
include { GATK4_GETPILEUPSUMMARIES        } from '../../../modules/nf-core/gatk4/getpileupsummaries'
include { GATK4_LEARNREADORIENTATIONMODEL } from '../../../modules/nf-core/gatk4/learnreadorientationmodel'
include { GATK4_MUTECT2                   } from '../../../modules/nf-core/gatk4/mutect2'

workflow BAM_SOMATIC_VARIANT_CALLING_GATK {
    take:
    ch_input // channel: [ val(meta), path(tumor_input), path(tumor_input_index), path(normal_input), path(normal_input_index)]
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

    ch_mutect2_in = ch_input.map{ meta, tumor_input, tumor_input_index, normal_input, normal_input_index ->
        if (normal_input && normal_input_index) {
            [meta, [tumor_input, normal_input], [tumor_input_index, normal_input_index]]
        }
        else {
            [meta, tumor_input, tumor_input_index]
        }
        }
        .combine(ch_interval_file)

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

    ch_mutect2_vcf = GATK4_MUTECT2.out.vcf
    ch_mutect2_stats = GATK4_MUTECT2.out.stats
    ch_mutect2_tbi = GATK4_MUTECT2.out.tbi
    ch_mutect2_f1r2 = GATK4_MUTECT2.out.f1r2

    // Generate artifactpriors using learnreadorientationmodel on the f1r2 output of GATK4_MUTECT2
    GATK4_LEARNREADORIENTATIONMODEL(ch_mutect2_f1r2)

    ch_artifact_priors = GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior

    // Generate pileup summary tables using getepileupsummaries
    ch_getpileupsummaries_in = ch_input
        .flatMap { meta, tumor_input, tumor_input_index, normal_input, normal_input_index ->
            def result = [
                [[id: "${meta.id}_tumor", original_meta: meta], // Add unique key so output files receive differing names
                tumor_input,
                tumor_input_index]
                ]

            if (normal_input && normal_input_index) {
                result << [
                    [id: "${meta.id}_normal", original_meta: meta], // Add unique key so output files receive differing names
                    normal_input,
                    normal_input_index
                    ]
            }

            return result
        }
        .combine(ch_interval_file)


    GATK4_GETPILEUPSUMMARIES(
        ch_getpileupsummaries_in,
        ch_fasta.collect(),
        ch_fai.map { meta, fai, _gzi -> [meta, fai] }.collect(),
        ch_dict.collect(),
        ch_germline_resource.collect(),
        ch_germline_resource_tbi.collect(),
    )

    ch_pileup_table = GATK4_GETPILEUPSUMMARIES.out.table
        .map { meta, table -> [meta.original_meta, table] } // Remove unique key and restore original meta

    // Contamination and segmentation tables created using calculatecontamination on the pileup summary table
    ch_calculatecontamination_in = ch_pileup_table
        .groupTuple()
        .map { meta, table_files ->
            [meta] + (table_files.size() > 1 ? table_files : [table_files[0], []])
            }

    GATK4_CALCULATECONTAMINATION(ch_calculatecontamination_in)

    ch_contamination_table = GATK4_CALCULATECONTAMINATION.out.contamination
    ch_segmentation_table = GATK4_CALCULATECONTAMINATION.out.segmentation

    // GATK4_MUTECT2 calls filtered by filtermutectcalls using the artifactpriors, contamination and segmentation tables
    ch_filtermutect_in = ch_mutect2_vcf
        .join(ch_mutect2_tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_mutect2_stats, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_artifact_priors, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_segmentation_table, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_contamination_table, failOnDuplicate: true, failOnMismatch: true)
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
    artifact_priors     = ch_artifact_priors // channel: [ val(meta), path(artifactprior) ]
    contamination_table = ch_contamination_table // channel: [ val(meta), path(table) ]
    filtered_stats      = GATK4_FILTERMUTECTCALLS.out.stats // channel: [ val(meta), path(stats) ]
    filtered_tbi        = GATK4_FILTERMUTECTCALLS.out.tbi // channel: [ val(meta), path(tbi) ]
    filtered_vcf        = GATK4_FILTERMUTECTCALLS.out.vcf // channel: [ val(meta), path(vcf) ]
    mutect2_f1r2        = ch_mutect2_f1r2 // channel: [ val(meta), path(f1r2) ]
    mutect2_stats       = ch_mutect2_stats // channel: [ val(meta), path(stats) ]
    mutect2_tbi         = ch_mutect2_tbi // channel: [ val(meta), path(tbi) ]
    mutect2_vcf         = ch_mutect2_vcf // channel: [ val(meta), path(vcf) ]
    pileup_table        = ch_pileup_table // channel: [ val(meta), path(table) ]
    segmentation_table  = ch_segmentation_table // channel: [ val(meta), path(table) ]
}
