include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow XENIUMRANGER_IMPORT_SEGMENTATION_IMAGE {

    take:
    ch_bundle                  // [meta, bundle]
    ch_nuclei                  // [meta, nuclei] (use [] if absent)
    ch_cells                   // [meta, cells]  (use [] if absent)
    ch_coordinate_transform    // [meta, coordinate_transform] (use [] if absent)

    main:

    def input = ch_bundle
        .map { meta, bundle -> [meta.id, [meta, bundle]] }
        .join(
            ch_nuclei.map { meta, nuclei -> [meta.id, [meta, nuclei]] }
            , by: 0, remainder: true)
        .join(
            ch_cells.map { meta, cells -> [meta.id, [meta, cells]] }
            , by: 0, remainder: true)
        .join(
            ch_coordinate_transform.map { meta, transform -> [meta.id, [meta, transform]] }
            , by: 0, remainder: true)
        .map { _id, bundle_entry, nuclei_entry, cells_entry, transform_entry ->

            def meta   = bundle_entry[0]
            def bundle = bundle_entry[1]

            def meta_n = nuclei_entry ? nuclei_entry[0] : [:]
            def nuclei = nuclei_entry ? nuclei_entry[1] : []

            def meta_c = cells_entry ? cells_entry[0] : [:]
            def cells  = cells_entry ? cells_entry[1] : []

            def meta_tr = transform_entry ? transform_entry[0] : [:]
            def transform = transform_entry ? transform_entry[1] : []

            def final_meta = meta_tr + meta_c + meta_n + meta

            [
                final_meta,
                bundle,
                [],
                [],
                nuclei,
                cells,
                transform
            ]
        }

    XENIUMRANGER_IMPORT_SEGMENTATION(input)

    emit:
    outs     = XENIUMRANGER_IMPORT_SEGMENTATION.out.outs
    versions = XENIUMRANGER_IMPORT_SEGMENTATION.out.versions_xeniumranger
}
