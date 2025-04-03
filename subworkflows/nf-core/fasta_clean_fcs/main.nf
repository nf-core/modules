include { FCS_FCSADAPTOR  } from '../../../modules/nf-core/fcs/fcsadaptor/main'
include { FCS_FCSGX       } from '../../../modules/nf-core/fcs/fcsgx/main'

workflow FASTA_CLEAN_FCS {

    take:
    ch_fasta // channel: [ val(meta), [ fasta ] ]
    database // value: the db to use

    main:
    ch_versions = Channel.empty()

    FCS_FCSADAPTOR ( ch_fasta )
    ch_versions = ch_versions.mix(FCS_FCSADAPTOR.out.versions)

    FCS_FCSGX ( FCS_FCSADAPTOR.out.cleaned_assembly, database )
    ch_versions = ch_versions.mix(FCS_FCSGX.out.versions)

    emit:
    fcsadaptor_cleaned_assembly = FCS_FCSADAPTOR.out.cleaned_assembly       // channel: [ val(meta), [ cleaned_assembly ] ]
    fcsadaptor_report           = FCS_FCSADAPTOR.out.adaptor_report         // channel: [ val(meta), [ adaptor_report ] ]
    fcsadaptor_log              = FCS_FCSADAPTOR.out.log                    // channel: [ val(meta), [ log ] ]
    fcsadaptor_pipeline_args    = FCS_FCSADAPTOR.out.pipeline_args          // channel: [ val(meta), [ pipeline_args ] ]
    fcsadaptor_skipped_trims    = FCS_FCSADAPTOR.out.skipped_trims          // channel: [ val(meta), [ skipped_trims ] ]
    fcsgx_report                = FCS_FCSGX.out.fcs_gx_report               // channel: [ val(meta), [ fcs_gx_report ] ]
    fcsgx_taxonomy_report       = FCS_FCSGX.out.taxonomy_report             // channel: [ val(meta), [ taxonomy_report ] ]
    versions                    = ch_versions                               // channel: [ versions.yml ]
}
