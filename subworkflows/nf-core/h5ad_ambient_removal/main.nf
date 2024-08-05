include { CELDA_DECONTX               } from '../../../../modules/nf-core/celda/decontx'
include { CELLBENDER_REMOVEBACKGROUND } from '../../../../modules/nf-core/cellbender/removebackground'
include { CELLBENDER_MERGE            } from '../../../../modules/nf-core/cellbender/merge'
include { SOUPX                       } from '../../../../modules/nf-core/soupx'
include { SCVITOOLS_SCAR              } from '../../../../modules/nf-core/scvitools/scar'

workflow H5AD_AMBIENT_REMOVAL {
    take:
    ch_filtered_unfiltered
    tool

    main:
    ch_versions = Channel.empty()

    if (tool == 'none') {
        println "AMBIENT_RNA_REMOVAL: Not performed since 'none' selected."
        ch_h5ad = ch_filtered_unfiltered.map{ meta, filtered, unfiltered -> [meta, filtered] }
    }
    else if (tool == 'decontx') {
        CELDA_DECONTX(ch_filtered_unfiltered)
        ch_h5ad = CELDA_DECONTX.out.h5ad
        ch_versions = ch_versions.mix(CELDA_DECONTX.out.versions)
    }
    else if (tool == 'cellbender') {
        CELLBENDER_REMOVEBACKGROUND(ch_filtered_unfiltered.map{meta, filtered, unfiltered -> [meta, unfiltered]})
        ch_versions = ch_versions.mix(CELLBENDER_REMOVEBACKGROUND.out.versions)

        CELLBENDER_MERGE(ch_filtered_unfiltered.map{ meta, filtered, raw -> [meta.id, meta, filtered, raw] }
            .join(CELLBENDER_REMOVEBACKGROUND.out.h5.map{ meta, h5 -> [meta.id, h5] }, by: 0, failOnMismatch: true)
            .map{ id, meta, filtered, raw, h5 -> [meta, filtered, raw, h5] })
        ch_h5ad = CELLBENDER_MERGE.out.h5ad
        ch_versions = ch_versions.mix(CELLBENDER_MERGE.out.versions)
    }
    else if (tool == 'soupx') {
        SOUPX(ch_filtered_unfiltered)
        ch_h5ad = SOUPX.out.h5ad
        ch_versions = ch_versions.mix(SOUPX.out.versions)
    }
    else if (tool == 'scar') {
        SCVITOOLS_SCAR(ch_filtered_unfiltered)
        ch_h5ad = SCVITOOLS_SCAR.out.h5ad
        ch_versions = SCVITOOLS_SCAR.out.versions
    }
    else {
        error "AMBIENT_REMOVAL: Unexpected tool selected: '${tool}'."
    }

    emit:
    h5ad = ch_h5ad

    versions = ch_versions
}
