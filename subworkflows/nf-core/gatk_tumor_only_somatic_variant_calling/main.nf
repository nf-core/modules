//
// Run GATK mutect2 in tumor only mode, getepileupsummaries, calculatecontamination and filtermutectcalls
//

params.mutect2_options     = [:]
params.getpileup_options   = [:]
params.calccontam_options  = [:]
params.filtercalls_options = [suffix: '_filtered']

include { GATK4_MUTECT2                as MUTECT2 }                  from '../../../modules/nf-core/gatk4/mutect2/main'
include { GATK4_GETPILEUPSUMMARIES     as GETPILEUPSUMMARIES }       from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION as CALCULATECONTAMINATION }   from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS      as FILTERMUTECTCALLS }        from '../../../modules/nf-core/gatk4/filtermutectcalls/main'

workflow GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING {
    take:
    input                     // channel: [ val(meta), path(input), path(input_index) ]
    fasta                     // channel: [ path(fasta)]                     /path/to/reference/fasta
    fai                       // channel: [ path(fai)]                       /path/to/reference/fasta/index
    dict                      // channel: [ path(dict)]                      /path/to/reference/fasta/dictionary
    germline_resource         // channel: [ path(resource)]                  /path/to/germline/resource
    germline_resource_tbi     // channel: [ path(resource_tbi)]              /path/to/germline/index
    panel_of_normals          // channel: [ path(pon)]                       /path/to/panel/of/normals
    panel_of_normals_tbi      // channel: [ path(index)]                     /path/to/panel/of/normals/index
    interval_file             // channel: [ path(interval)]                  /path/to/interval/file


    main:
    ch_versions = Channel.empty()
    mutect2_input = channel.from(input)

    //
    //Perform variant calling using mutect2 module in tumor single mode.
    //
    MUTECT2 ( mutect2_input , true , false , false , [] , fasta , fai , dict , germline_resource , germline_resource_tbi , panel_of_normals , panel_of_normals_tbi )
    ch_versions = ch_versions.mix(MUTECT2.out.versions)

    //
    //Generate pileup summary table using getepileupsummaries.
    //
    pileup_input = channel.from(input).map {
        meta, input_file, input_index, which_norm ->
        [meta, input_file[0], input_index[0]]
    }
    GETPILEUPSUMMARIES ( pileup_input , germline_resource , germline_resource_tbi , interval_file )
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES.out.versions)

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
    ch_filtermutect_in = ch_vcf.combine(ch_tbi, by: 0).combine(ch_stats, by: 0).combine(ch_segment, by: 0).combine(ch_contamination, by: 0)
    FILTERMUTECTCALLS ( ch_filtermutect_in, fasta, fai, dict )
    ch_versions = ch_versions.mix(FILTERMUTECTCALLS.out.versions)

    emit:
    mutect2_vcf            = MUTECT2.out.vcf.collect()                             // channel: [ val(meta), path(vcf)]
    mutect2_index          = MUTECT2.out.tbi.collect()                             // channel: [ val(meta), path(tbi)]
    mutect2_stats          = MUTECT2.out.stats.collect()                           // channel: [ val(meta), path(stats)]

    pileup_table           = GETPILEUPSUMMARIES.out.table.collect()                // channel: [ val(meta), path(table)]

    contamination_table    = CALCULATECONTAMINATION.out.contamination.collect()    // channel: [ val(meta), path(contamination) ]
    segmentation_table     = CALCULATECONTAMINATION.out.segmentation.collect()     // channel: [ val(meta), path(segmentation) ]

    filtered_vcf           = FILTERMUTECTCALLS.out.vcf.collect()                   // channel: [ val(meta), path(vcf) ]
    filtered_index         = FILTERMUTECTCALLS.out.tbi.collect()                   // channel: [ val(meta), path(tbi) ]
    filtered_stats         = FILTERMUTECTCALLS.out.stats.collect()                 // channel: [ val(meta), path(stats)]

    versions               = ch_versions                                           // channel: [ path(versions.yml) ]
}
