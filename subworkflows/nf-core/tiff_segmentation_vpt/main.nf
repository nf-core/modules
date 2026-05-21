//
// Perform segmentation of MERSCOPE TIFF images using the vizgen-postprocessing tool
//

include { VIZGENPOSTPROCESSING_PREPARESEGMENTATION     } from '../../../modules/nf-core/vizgenpostprocessing/preparesegmentation/main'
include { VIZGENPOSTPROCESSING_RUNSEGMENTATIONONTILE   } from '../../../modules/nf-core/vizgenpostprocessing/runsegmentationontile/main'
include { VIZGENPOSTPROCESSING_COMPILETILESEGMENTATION } from '../../../modules/nf-core/vizgenpostprocessing/compiletilesegmentation/main'

workflow TIFF_SEGMENTATION_VPT {

    take:
    ch_input             // channel: [ meta, image_dir, micron_to_mosaic_pixel_transform.csv ]
    ch_algorithm_json    // channel: [ algorithm.json ]
    ch_images_regex      // channel: [ val(regex) ]
    ch_custom_weights    // channel: [ custom_weights ]

    main:

    ch_versions = channel.empty()

    //
    // Generate segmentation parameters file
    //
    VIZGENPOSTPROCESSING_PREPARESEGMENTATION (
        ch_input,
        ch_algorithm_json,
        ch_images_regex
    )
    ch_versions = ch_versions.mix(VIZGENPOSTPROCESSING_PREPARESEGMENTATION.out.versions)

    // Create list sequence of 0..N tiles
    // Get N from segmentation params JSON file
    VIZGENPOSTPROCESSING_PREPARESEGMENTATION.out.segmentation_files
        .map { _meta, seg_params -> seg_params }
        .splitJson(path: 'window_grid' )
        .filter { index -> index.key == 'num_tiles' }
        .map { index -> index.value.toInteger() }
        .flatMap { num -> (0..num-1).toList() }
        .set{ ch_tiles }

    // Create channel containing required segmentation files
    VIZGENPOSTPROCESSING_PREPARESEGMENTATION.out.segmentation_files
        .join(ch_input, by: 0)
        .map { meta, seg_params, image_dir, _pixel_file ->
            tuple(meta, image_dir, seg_params)
        }
        .set { ch_segmentation_files }

    //
    // Run vpt segmentation on each tile
    //
    VIZGENPOSTPROCESSING_RUNSEGMENTATIONONTILE(
        ch_segmentation_files.combine(ch_tiles),
        ch_algorithm_json,
        ch_custom_weights
    )
    ch_versions = ch_versions.mix(VIZGENPOSTPROCESSING_RUNSEGMENTATIONONTILE.out.versions.first())

    /// collect segmented tiles
    VIZGENPOSTPROCESSING_RUNSEGMENTATIONONTILE.out.segmented_tile
        .map { _meta, seg_tile -> seg_tile }
        .flatten()
        .collect()
        .set{ ch_segmented_tiles }

    //
    // Compile segmentation results from all tiles
    //
    VIZGENPOSTPROCESSING_COMPILETILESEGMENTATION(
        ch_segmentation_files,
        ch_algorithm_json,
        ch_segmented_tiles
    )
    ch_versions = ch_versions.mix(VIZGENPOSTPROCESSING_COMPILETILESEGMENTATION.out.versions.first())

    // Output channels for segmentations
    ch_micron_space_parquet =
        VIZGENPOSTPROCESSING_COMPILETILESEGMENTATION.out.micron_space
    ch_mosaic_space_parquet =
        VIZGENPOSTPROCESSING_COMPILETILESEGMENTATION.out.mosaic_space

    emit:
    micron_space_segmentation       = ch_micron_space_parquet  // channel: [ meta, parquet ]
    mosaic_space_segmentation       = ch_mosaic_space_parquet  // channel: [ meta, parquet ]
    versions                        = ch_versions              // channel: [ versions.yml ]
}
