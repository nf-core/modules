include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow XENIUMRANGER_IMPORT_SEGMENTATION_TRANSCRIPT {

    take:
    ch_bundle                // [meta, bundle]
    ch_transcript_assignment // [meta, transcript_assignment]
    ch_viz_polygons          // [meta, viz_polygons]

    main:

    def input = ch_bundle
        .map { meta, bundle -> [meta.id, meta, bundle] }
        .combine(
            ch_transcript_assignment
                .map { meta, transcript_assignment -> [meta.id, meta, transcript_assignment] }
            , by:0)
        .combine(
            ch_viz_polygons
                .map { meta, viz_polygons -> [meta.id, meta, viz_polygons] }
            , by:0)
        .map { _id, meta, bundle, meta_t, transcript_assignment, meta_v, viz_polygons -> [
            meta_t + meta_v + meta,
            bundle,
            transcript_assignment,
            viz_polygons,
            [],
            [],
            []
            ]
        }

    XENIUMRANGER_IMPORT_SEGMENTATION(input)

    emit:
    outs     = XENIUMRANGER_IMPORT_SEGMENTATION.out.outs
    versions = XENIUMRANGER_IMPORT_SEGMENTATION.out.versions_xeniumranger
}
