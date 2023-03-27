//
// Run GATK mutect2 in tumor only mode, getepileupsummaries, calculatecontamination and filtermutectcalls
//

params.mutect2_options     = [:]
params.getpileup_options   = [:]
params.calccontam_options  = [:]
params.filtercalls_options = [suffix: '_filtered']

include { GATK4_MUTECT2                   as MUTECT2         } from '../../../modules/nf-core/gatk4/mutect2/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES        } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION as CALCULATECONTAMINATION }   from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS      as FILTERMUTECTCALLS }        from '../../../modules/nf-core/gatk4/filtermutectcalls/main'

workflow BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK  {
    take:

    ch_input                     // channel: [ val(meta), path(input), path(input_index), [] ]
    ch_fasta                     // channel: [ path(fasta) ]
    ch_fai                       // channel: [ path(fai) ]
    ch_dict                      // channel: [ path(dict) ]
    ch_germline_resource         // channel: [ path(germline_resource) ]
    ch_germline_resource_tbi     // channel: [ path(germline_resource_tbi) ]
    ch_panel_of_normals          // channel: [ path(panel_of_normals) ]
    ch_panel_of_normals_tbi      // channel: [ path(panel_of_normals_tbi) ]
    ch_interval_file             // channel: [ path(interval_file) ]


    main:
    ch_versions = Channel.empty()
    germline_resource_pileup     = ch_germline_resource_tbi ? germline_resource : Channel.empty()
    germline_resource_pileup_tbi = ch_germline_resource_tbi ?: Channel.empty()


    //
    //Perform variant calling using mutect2 module in tumor single mode.
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
    //Generate pileup summary table using getepileupsummaries.
    //
    ch_pileup_tumor_input = ch_input.map {
        meta, input_file, input_index, which_norm ->
        [meta, input_file[0], input_index[0]]
    }

    GETPILEUPSUMMARIES(ch_pileup_tumor_input,
            fasta,
            fai,
            dict,
            germline_resource_pileup,
            germline_resource_pileup_tbi)

    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES.out.versions.first())

    //
    //Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //
    ch_pileup = GETPILEUPSUMMARIES.out.table.collect()
    //[] is a placeholder for the optional input where the matched normal sample would be passed in for tumor-normal samples, which is not necessary for this workflow.
    ch_pileup.add([])
    CALCULATECONTAMINATION ( ch_pileup, true )
    ch_versions = ch_versions.mix(CALCULATECONTAMINATION.out.versions)

    //
    //Mutect2 calls filtered by filtermutectcalls using the contamination and segmentation tables.
    //
    ch_vcf =           MUTECT2.out.vcf.collect()
    ch_tbi =           MUTECT2.out.tbi.collect()
    ch_stats =         MUTECT2.out.stats.collect()
    //[] is added as a placeholder for the optional input file artifact priors, which is only used for tumor-normal samples and therefor isn't needed in this workflow.
    ch_stats.add([])
    ch_segment =       CALCULATECONTAMINATION.out.segmentation.collect()
    ch_contamination = CALCULATECONTAMINATION.out.contamination.collect()
    //[] is added as a placeholder for entering a contamination estimate value, which is not needed as this workflow uses the contamination table instead.
    ch_contamination.add([])

    ch_filtermutect_in = ch_vcf
        .join(ch_tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_stats, failOnDuplicate: true, failOnMismatch: true)
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

    pileup_table           = GETPILEUPSUMMARIES.out.table.collect()                // channel: [ val(meta), path(table) ]

    contamination_table    = CALCULATECONTAMINATION.out.contamination.collect()    // channel: [ val(meta), path(table) ]
    segmentation_table     = CALCULATECONTAMINATION.out.segmentation.collect()     // channel: [ val(meta), path(table) ]

    filtered_vcf           = FILTERMUTECTCALLS.out.vcf.collect()                   // channel: [ val(meta), path(vcf) ]
    filtered_tbi           = FILTERMUTECTCALLS.out.tbi.collect()                   // channel: [ val(meta), path(tbi) ]
    filtered_stats         = FILTERMUTECTCALLS.out.stats.collect()                 // channel: [ val(meta), path(stats) ]

    versions               = ch_versions                                           // channel: [ path(versions.yml) ]
}
