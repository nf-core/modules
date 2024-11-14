//
// cellbender remove background and empty droplet detection and removal
//
include { CELLBENDER_REMOVEBACKGROUND } from '../../../modules/nf-core/cellbender/removebackground'
include { ANNDATA_BARCODES            } from '../../../modules/nf-core/anndata/barcodes'

workflow CELLBENDER_EMPTY_DROPLET_REMOVAL {

    take:
    ch_unfiltered

    main:
    ch_versions = Channel.empty()

    CELLBENDER_REMOVEBACKGROUND(ch_unfiltered)
    ch_versions = ch_versions.mix(CELLBENDER_REMOVEBACKGROUND.out.versions)

    ch_combined = ch_unfiltered.join(CELLBENDER_REMOVEBACKGROUND.out.barcodes)

    ANNDATA_BARCODES(ch_combined)
    ch_versions = ch_versions.mix(ANNDATA_BARCODES.out.versions)

    ch_h5ad = ANNDATA_BARCODES.out.h5ad

    emit:
    h5ad = ch_h5ad

    versions = ch_versions
}

