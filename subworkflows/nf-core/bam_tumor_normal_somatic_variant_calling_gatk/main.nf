//
// Run GATK mutect2 in tumor normal mode, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls
//

include { GATK4_MUTECT2                   as MUTECT2 }                   from '../../../modules/nf-core/gatk4/mutect2/main'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNREADORIENTATIONMODEL } from '../../../modules/nf-core/gatk4/learnreadorientationmodel/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_TUMOR }  from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_NORMAL}  from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION    as CALCULATECONTAMINATION }    from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS         as FILTERMUTECTCALLS }         from '../../../modules/nf-core/gatk4/filtermutectcalls/main'

workflow BAM_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING_GATK {
    take:
    ch_input                 // channel: [ val(meta), path(input), path(input_index), val(which_norm) ]
    ch_fasta                 // channel: /path/to/reference/fasta
    ch_fai                   // channel: /path/to/reference/fasta/index
    ch_dict                  // channel: /path/to/reference/fasta/dictionary
    ch_germline_resource     // channel: /path/to/germline/resource
    ch_germline_resource_tbi // channel: /path/to/germline/index
    ch_panel_of_normals      // channel: /path/to/panel/of/normals
    ch_panel_of_normals_tbi  // channel: /path/to/panel/of/normals/index
    ch_interval_file         // channel: /path/to/interval/file

    main:
    ch_versions = Channel.empty()

    //
    // Perform variant calling using mutect2 module in tumor single mode.
    //
    MUTECT2 (
        ch_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_tbi,
        ch_panel_of_normals,
        ch_panel_of_normals_tbi
    )

    ch_versions = ch_versions.mix(MUTECT2.out.versions)

    //
    // Generate artifactpriors using learnreadorientationmodel on the f1r2 output of mutect2.
    //
    LEARNREADORIENTATIONMODEL (MUTECT2.out.f1r2.collect())
    ch_versions = ch_versions.mix(LEARNREADORIENTATIONMODEL.out.versions)

    //
    // Generate pileup summary tables using getepileupsummaries. tumor sample should always be passed in as the first input and input list entries of ch_mutect2_in,
    // to ensure correct file order for calculatecontamination.
    //
    ch_pileup_tumor_input = ch_input.combine(ch_interval_file).map {
        meta, input_file, input_index, which_norm, intervals ->
        [meta, input_file[0], input_index[0], intervals]
    }

    ch_pileup_normal_input = ch_input.combine(ch_interval_file).map {
        meta, input_file, input_index, which_norm, intervals ->
        [meta, input_file[1], input_index[1], intervals]
    }

    GETPILEUPSUMMARIES_TUMOR (
        ch_pileup_tumor_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_tbi
    )

    GETPILEUPSUMMARIES_NORMAL (
        ch_pileup_normal_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_tbi
    )

    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_TUMOR.out.versions.first())
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_NORMAL.out.versions.first())

    //
    // Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //
    ch_pileup_tumor = GETPILEUPSUMMARIES_TUMOR.out.table.collect()
    ch_pileup_normal = GETPILEUPSUMMARIES_NORMAL.out.table.collect()
    ch_calccon_in = ch_pileup_tumor.join(ch_pileup_normal, failOnDuplicate: true, failOnMismatch: true)
    CALCULATECONTAMINATION ( ch_calccon_in )
    ch_versions   = ch_versions.mix(CALCULATECONTAMINATION.out.versions)

    //
    // Mutect2 calls filtered by filtermutectcalls using the artifactpriors, contamination and segmentation tables.
    //
    ch_vcf           = MUTECT2.out.vcf.collect()
    ch_tbi           = MUTECT2.out.tbi.collect()
    ch_stats         = MUTECT2.out.stats.collect()
    ch_orientation   = LEARNREADORIENTATIONMODEL.out.artifactprior.collect()
    ch_segment       = CALCULATECONTAMINATION.out.segmentation.collect()
    ch_contamination = CALCULATECONTAMINATION.out.contamination.collect()

    //[] is used as a placeholder for optional input to specify the contamination estimate as a value, since the contamination table is used, this is not needed.
    ch_contamination.add([])
    ch_filtermutect_in = ch_vcf
        .join(ch_tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_stats, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_orientation, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_segment, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_contamination, failOnDuplicate: true, failOnMismatch: true)

    FILTERMUTECTCALLS (
        ch_filtermutect_in,
        ch_fasta,
        ch_fai,
        ch_dict
    )
    ch_versions = ch_versions.mix(FILTERMUTECTCALLS.out.versions.first())

    emit:
    mutect2_vcf            = MUTECT2.out.vcf.collect()                             // channel: [ val(meta), path(vcf) ]
    mutect2_tbi            = MUTECT2.out.tbi.collect()                             // channel: [ val(meta), path(tbi) ]
    mutect2_stats          = MUTECT2.out.stats.collect()                           // channel: [ val(meta), path(stats) ]
    mutect2_f1r2           = MUTECT2.out.f1r2.collect()                            // channel: [ val(meta), path(f1r2) ]

    artifact_priors        = LEARNREADORIENTATIONMODEL.out.artifactprior.collect() // channel: [ val(meta), path(artifactprior) ]

    pileup_table_tumor     = GETPILEUPSUMMARIES_TUMOR.out.table.collect()          // channel: [ val(meta), path(table) ]
    pileup_table_normal    = GETPILEUPSUMMARIES_NORMAL.out.table.collect()         // channel: [ val(meta), path(table) ]

    contamination_table    = CALCULATECONTAMINATION.out.contamination.collect()    // channel: [ val(meta), path(table) ]
    segmentation_table     = CALCULATECONTAMINATION.out.segmentation.collect()     // channel: [ val(meta), path(table) ]

    filtered_vcf           = FILTERMUTECTCALLS.out.vcf.collect()                   // channel: [ val(meta), path(vcf) ]
    filtered_tbi           = FILTERMUTECTCALLS.out.tbi.collect()                   // channel: [ val(meta), path(tbi) ]
    filtered_stats         = FILTERMUTECTCALLS.out.stats.collect()                 // channel: [ val(meta), path(stats) ]

    versions               = ch_versions                                           // channel: [ path(versions.yml) ]
}
