process RAMI2D_TRANSFORM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/009a819a706893896160b3ea33fcf25254d705f7bdccedf2baa63bcdd009c826/data':
        'community.wave.seqera.io/library/rami2d-env:80a2c0ef1e678902' }"

    input:
    tuple val(meta), path(image_file), path(transform_dir)
    tuple val(meta2), path(markers)
    val(mpp)

    output:
    tuple val(meta), path("${prefix}/*transformed.ome.tif"), emit: transformed_image
    tuple val("${task.process}"), val('rami2d'), eval('rami2d-transform -v | sed "s/rami2d-transform //"'), emit: versions_rami2d, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def markers_arg = markers ? "-m ${markers}" : ""

    """
    export MPLCONFIGDIR=./matplotlib_cache
    export XDG_CACHE_HOME=./.cache
    export FONTCONFIG_PATH=./.cache/fontconfig

    mkdir -p $prefix

    rami2d-transform \\
        -i $image_file \\
        -mpp $mpp \\
        -tdir $transform_dir \\
        -o $prefix \\
        $markers_arg \\
        $args
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p $prefix
    touch ${prefix}/${image_file.baseName}_transformed.ome.tif
    """
}
