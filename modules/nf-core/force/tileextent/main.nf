// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FORCE_TILEEXTENT {
    tag "${aoi.simpleName}"
    label 'process_single'
    stageInMode 'copy'
    // needed by the module to work properly when aoi is a shapefile

    container "quay.io/nf-core/force:3.8.01"

    input:
    path aoi
    path datacube_definition, stageAs: "tmp/datacube-definition.prj"
    path shapefile_dbf
    path shapefile_prj
    path shapefile_shx

    output:
    path "tile_allow.txt", emit: tile_allow
    tuple val("${task.process}"), val('force'), eval('force-tile-extent -v'), emit: versions_force, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    force-tile-extent -d tmp/ -a tile_allow.txt ${aoi}
    rm -r tmp
    """

    stub:
    """
    touch tile_allow.txt
    """
}
