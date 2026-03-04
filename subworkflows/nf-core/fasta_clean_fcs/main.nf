include { FCS_FCSADAPTOR  } from '../../../modules/nf-core/fcs/fcsadaptor/main'
include { FCSGX_RUNGX     } from '../../../modules/nf-core/fcsgx/rungx/main'

workflow FASTA_CLEAN_FCS {

    take:
    ch_fasta     // channel: [ val(meta), [ fasta ] ]
    database     // channel: [ db ]
    ramdisk_path // value: the path to the ramdisk

    main:
    ch_versions = channel.empty()

    ch_fasta.map{
        meta, _fasta -> [ meta.taxid ?: error("taxid is mandatory in the meta map") ]
    }

    FCS_FCSADAPTOR ( ch_fasta )
    ch_versions = ch_versions.mix(FCS_FCSADAPTOR.out.versions)

    ch_cleaned_assembly = ch_fasta
        .join(FCS_FCSADAPTOR.out.cleaned_assembly, by:0, remainder: true )
        .map { meta, input, cleaned ->
            [ meta, meta.taxid, cleaned ?: input ]
        }

    FCSGX_RUNGX ( ch_cleaned_assembly, database, ramdisk_path)
    ch_versions = ch_versions.mix(FCSGX_RUNGX.out.versions)

    emit:
    fcsadaptor_cleaned_assembly = ch_cleaned_assembly                       // channel: [ val(meta), [ cleaned_assembly ] ]
    fcsadaptor_report           = FCS_FCSADAPTOR.out.adaptor_report         // channel: [ val(meta), [ adaptor_report ] ]
    fcsadaptor_log              = FCS_FCSADAPTOR.out.log                    // channel: [ val(meta), [ log ] ]
    fcsadaptor_pipeline_args    = FCS_FCSADAPTOR.out.pipeline_args          // channel: [ val(meta), [ pipeline_args ] ]
    fcsadaptor_skipped_trims    = FCS_FCSADAPTOR.out.skipped_trims          // channel: [ val(meta), [ skipped_trims ] ]
    fcsgx_report                = FCSGX_RUNGX.out.fcsgx_report              // channel: [ val(meta), [ fcs_gx_report ] ]
    fcsgx_taxonomy_report       = FCSGX_RUNGX.out.taxonomy_report           // channel: [ val(meta), [ taxonomy_report ] ]
    fcsgx_log                   = FCSGX_RUNGX.out.log                       // channel: [ val(meta), [ log ] ]
    fcsgx_hits                  = FCSGX_RUNGX.out.hits                      // channel: [ val(meta), [ hits ] ]
    versions                    = ch_versions                               // channel: [ versions.yml ]
}
