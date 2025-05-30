#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process PREPARE_TEST_DATA {
    input:
    val dummy_input

    output:
    path "images", emit: image_directory
    path "images/micron_to_mosaic_pixel_transform.csv", emit: micron_to_mosaic_transform
    path "detected_transcripts.csv", emit: transcripts_csv

    when:

    script:
    """
    # Set up images and directory structure
    mkdir -p images
    wget -O images/mosaic_DAPI_z3.tif https://github.com/nf-core/test-datasets/raw/modules/data/imaging/segmentation/nuclear_image.tif
    cp images/mosaic_DAPI_z3.tif images/mosaic_PolyT_z3.tif

    # Create a dummy micron to mosaic pixel transform file
    echo "2 0 0\n0 2 0\n0 0 1" > images/micron_to_mosaic_pixel_transform.csv

    # Create a dummy detected transcripts CSV file
    echo ",barcode_id,global_x,global_y,global_z,x,y,fov,gene,transcript_id" > detected_transcripts.csv
    echo "101,10,370.95007,5.520504,0.0,611.833,110.285774,0,GAPDH,ENST00000542121.2" >> detected_transcripts.csv
    echo "102,10,355.46716,6.4616404,0.0,668.4725,119.0,0,GAPDH,ENST00000542121.2" >> detected_transcripts.csv
    """
}
