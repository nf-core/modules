process FORCE_TILEEXTENT {
    tag "${aoi.simpleName}"
    label 'process_single'
    stageInMode 'copy' // needed by the module to work properly when aoi is a shapefile

    container "nf-core/force:3.8.01"

    input:
    path aoi
    path(datacube_definition, stageAs: "tmp/datacube-definition.prj")
    path shapefile_dbf
    path shapefile_prj
    path shapefile_shx

    output:
    path "tile_allow.txt", emit: tile_allow
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    force-tile-extent -d tmp/ -a tile_allow.txt $aoi
    rm -r tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        force: \$(force-tile-extent -v)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch tile_allow.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        force: \$(force-tile-extent -v)
    END_VERSIONS
    """
}
