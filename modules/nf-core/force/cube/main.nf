process FORCE_CUBE {
    tag "${aoi.simpleName}"
    label 'process_single'

    container "nf-core/force:3.8.01"

    input:
    path aoi
    path 'mask/datacube-definition.prj'
    path shapefile_dbf
    path shapefile_prj
    path shapefile_shx

    output:
    path 'mask/*/*.tif', emit: masks
    tuple val("${task.process}"), val('force'), eval('force-cube -v'), emit: versions_force, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    force-cube $args -o mask/ $aoi
    """

    stub:
    def args     = task.ext.args ?: ''
    def matcher  = args =~ /-b\s+(\S+)/
    def baseName = args.contains('-b') && matcher ? matcher[0][1] : aoi.getSimpleName()

    def tile1 = "X0000_Y0000"
    def tile2 = "X0001_Y0000"
    """
    mkdir mask/$tile1 mask/$tile2
    touch mask/$tile1/${baseName}.tif mask/$tile2/${baseName}.tif
    """
}
