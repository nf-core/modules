//
// Register H&E stained and multiplexed tissue images and transform segmentation masks using stainwarpy
//

include { STAINWARPY_REGISTER         } from '../../../modules/nf-core/stainwarpy/register/main'
include { STAINWARPY_EXTRACTCHANNEL   } from '../../../modules/nf-core/stainwarpy/extractchannel/main'
include { STAINWARPY_TRANSFORMSEGMASK } from '../../../modules/nf-core/stainwarpy/transformsegmask/main'

workflow TIF_REGISTRATION_STAINWARPY {

    take:
    ch_hne              // channel: [ val(meta), path to .tif ]
    ch_multiplexed      // channel: [ val(meta), path to .tif ]
    ch_segmask          // channel: [ val(meta), path to .tif ] (optional)
    val_fixed_img       // val: fixed image to use ("multiplexed" or "hne")
    val_final_img_sz    // val: final image size to use ("multiplexed" or "hne")

    main:
    ch_transformed_segmask = channel.empty()
    ch_multiplexed_forward = channel.empty()

    if ( val_fixed_img == 'multiplexed') {
        STAINWARPY_EXTRACTCHANNEL ( ch_multiplexed )
        ch_multiplexed_forward = STAINWARPY_EXTRACTCHANNEL.out.single_ch_image
        STAINWARPY_REGISTER ( ch_hne, ch_multiplexed_forward,  val_fixed_img, val_final_img_sz )
    } else {
        STAINWARPY_REGISTER ( ch_hne, ch_multiplexed, val_fixed_img, val_final_img_sz )
        ch_multiplexed_forward = ch_multiplexed
    }

    if ( ch_segmask ) {
        STAINWARPY_TRANSFORMSEGMASK ( ch_hne, ch_multiplexed_forward, ch_segmask, STAINWARPY_REGISTER.out.tform_map, val_fixed_img, val_final_img_sz )
        ch_transformed_segmask = STAINWARPY_TRANSFORMSEGMASK.out.transformed_seg_mask
    }

    emit:
    transformed_image   = STAINWARPY_REGISTER.out.reg_image                     // channel: [ val(meta), *_transformed_image.ome.tif             ]
    metrics             = STAINWARPY_REGISTER.out.reg_metrics_tform             // channel: [ val(meta), *_registration_metrics_tform_map.json   ]
    transformed_segmask = ch_transformed_segmask                                // channel: [ val(meta), *_transformed_segmentation_mask.ome.tif ]
}
