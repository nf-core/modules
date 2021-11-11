//
// Run GATK mutect2 in tumour only mode, getepileupsummaries, calculatecontamination and filtermutectcalls
//

params.mutect2_options     = [:]
params.getpileup_options   = [:]
params.calccontam_options  = [:]
params.filtercalls_options = [suffix: '_filtered']

include { GATK4_MUTECT2                     } from '../../../modules/gatk4/mutect2/main'                addParams( options: params.mutect2_options )
include { GATK4_GETPILEUPSUMMARIES          } from '../../../modules/gatk4/getpileupsummaries/main'    addParams( options: params.getpileup_options )
include { GATK4_CALCULATECONTAMINATION      } from '../../../modules/gatk4/calculatecontamination/main' addParams( options: params.calccontam_options )
include { GATK4_FILTERMUTECTCALLS           } from '../../../modules/gatk4/filtermutectcalls/main'      addParams( options: params.filtercalls_options )

workflow GATK_TUMOUR_ONLY_SOMATIC_VARIANT_CALLING {
    take:
    ch_mutect2_in             // channel: [ val(meta), [ input ], [ input_index ], [] ]
    fasta                     // channel: /path/to/reference/fasta
    fastaidx                  // channel: /path/to/reference/fasta/index
    dict                      // channel: /path/to/reference/fasta/dictionary
    germline_resource         // channel: /path/to/germline/resource
    germline_resource_idx     // channel: /path/to/germline/index
    panel_of_normals          // channel: /path/to/panel/of/normals
    panel_of_normals_idx      // channel: /path/to/panel/of/normals/index
    interval_file             // channel: /path/to/interval/file


    main:
    ch_versions      = Channel.empty()
    input = channel.from(ch_mutect2_in)
    //
    //Perform variant calling using mutect2 module in tumour single mode.
    //
    GATK4_MUTECT2 ( input , true , false , false , [] , fasta , fastaidx , dict , germline_resource , germline_resource_idx , panel_of_normals , panel_of_normals_idx )
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions)

    //
    //Generate pileup summary table using getepileupsummaries.
    //
    pileup_input = channel.from(ch_mutect2_in)
    pileup_input = pileup_input.map {
        meta, input_file, input_index, which_norm ->
        [meta, input_file[0], input_index[0]]
    }
    GATK4_GETPILEUPSUMMARIES ( pileup_input , germline_resource , germline_resource_idx , interval_file )
    ch_versions = ch_versions.mix(GATK4_GETPILEUPSUMMARIES.out.versions)

    //
    //Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //
    ch_pileup = GATK4_GETPILEUPSUMMARIES.out.table.collect()
    ch_pileup.add([])
    GATK4_CALCULATECONTAMINATION ( ch_pileup, true )
    ch_versions = ch_versions.mix(GATK4_CALCULATECONTAMINATION.out.versions)

    //
    //Mutect2 calls filtered by filtermutectcalls using the contamination and segmentation tables.
    //
    ch_vcf = GATK4_MUTECT2.out.vcf.collect()
    ch_tbi = GATK4_MUTECT2.out.tbi.collect()
    ch_stats = GATK4_MUTECT2.out.stats.collect()
    ch_stats.add([])
    ch_segment = GATK4_CALCULATECONTAMINATION.out.segmentation.collect()
    ch_contamination = GATK4_CALCULATECONTAMINATION.out.contamination.collect()
    ch_contamination.add([])
    ch_filtermutect_in = ch_vcf.combine(ch_tbi, by: 0).combine(ch_stats, by: 0).combine(ch_segment, by: 0).combine(ch_contamination, by: 0)
    GATK4_FILTERMUTECTCALLS ( ch_filtermutect_in, fasta, fastaidx, dict )
    ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS.out.versions)

    emit:
    mutect2_vcf            = GATK4_MUTECT2.out.vcf.collect()                             // channel: [ val(meta), [ vcf ] ]
    mutect2_index          = GATK4_MUTECT2.out.tbi.collect()                             // channel: [ val(meta), [ tbi ] ]
    mutect2_stats          = GATK4_MUTECT2.out.stats.collect()                           // channel: [ val(meta), [ stats ] ]

    pileup_table           = GATK4_GETPILEUPSUMMARIES.out.table.collect()                // channel: [ val(meta), [ table ] ]

    contamination_table    = GATK4_CALCULATECONTAMINATION.out.contamination.collect()    // channel: [ val(meta), [ contamination ] ]
    segmentation_table     = GATK4_CALCULATECONTAMINATION.out.segmentation.collect()     // channel: [ val(meta), [ segmentation ] ]

    filtered_vcf           = GATK4_FILTERMUTECTCALLS.out.vcf.collect()                  // channel: [ val(meta), [ vcf ] ]
    filtered_index         = GATK4_FILTERMUTECTCALLS.out.tbi.collect()                   // channel: [ val(meta), [ tbi ] ]
    filtered_stats         = GATK4_FILTERMUTECTCALLS.out.stats.collect()                 // channel: [ val(meta), [ stats ] ]

    versions               = ch_versions                                                 // channel: [ versions.yml ]
}
