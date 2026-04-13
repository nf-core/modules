//
// Apply cellbender and anndata to h5ad for background and empty droplet removal
//
include { CELLBENDER_REMOVEBACKGROUND } from '../../../modules/nf-core/cellbender/removebackground'
include { ANNDATA_BARCODES            } from '../../../modules/nf-core/anndata/barcodes'

workflow H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA {

    take:
    ch_unfiltered // channel: [mandatory] meta, h5ad

    main:
    CELLBENDER_REMOVEBACKGROUND(ch_unfiltered)

    ANNDATA_BARCODES(ch_unfiltered.join(CELLBENDER_REMOVEBACKGROUND.out.barcodes))

    emit:
    h5ad = ANNDATA_BARCODES.out.h5ad  // channel: [ val(meta), path(h5ad) ]
}
