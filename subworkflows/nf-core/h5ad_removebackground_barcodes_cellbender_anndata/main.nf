//
// Apply cellbender and anndata to h5ad for background and empty droplet removal
//
include { CELLBENDER_REMOVEBACKGROUND } from '../../../modules/nf-core/cellbender/removebackground'
include { ANNDATA_BARCODES            } from '../../../modules/nf-core/anndata/barcodes'

workflow H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA {

    take:
    ch_unfiltered // channel: [mandatory] meta, h5ad

    main:
    ch_versions = Channel.empty()

    CELLBENDER_REMOVEBACKGROUND(ch_unfiltered)
    ch_versions = ch_versions.mix(CELLBENDER_REMOVEBACKGROUND.out.versions)

    ch_combined = ch_unfiltered.join(CELLBENDER_REMOVEBACKGROUND.out.barcodes)

    ANNDATA_BARCODES(ch_combined)
    ch_versions = ch_versions.mix(ANNDATA_BARCODES.out.versions)

    emit:
    h5ad = ANNDATA_BARCODES.out.h5ad  // channel: [ val(meta), path(h5ad) ]

    versions = ch_versions  // channel: [ path(versions.yml) ]
}

