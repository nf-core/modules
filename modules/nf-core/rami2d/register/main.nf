process RAMI2D_REGISTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cd/cde38bd19f190e920f5f99178fd6802f841e5e52f596e2efe13f95c733637850/data':
        'community.wave.seqera.io/library/pip_python_tifffile_ome-types_pruned:d3b5c8c723a6348d' }"

    input:
    tuple val(meta), path(fixed_img), path(moving_img)
    tuple val(meta2), path(markers)
    val(ifix)
    val(imov)
    val(mpp_fix)
    val(mpp_mov)
    val(mpp_reg)

    output:
    tuple val(meta), path("${prefix}/*.ome.tif"), emit: registered_image
    tuple val(meta), path("${prefix}/qc")       , emit: qc_dir
    tuple val("${task.process}"), val('rami2d') , eval('rami2d-register -v | sed "s/rami2d-register //"'), emit: versions_rami2d, topic: versions

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

    rami2d-register \\
        -fix $fixed_img \\
        -mov $moving_img \\
        -ifix $ifix \\
        -imov $imov \\
        -mpp-fix $mpp_fix \\
        -mpp-mov $mpp_mov \\
        -mpp-reg $mpp_reg \\
        -o $prefix \\
        $markers_arg \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p $prefix/qc/refchns
    mkdir -p $prefix/qc/keypoints
    mkdir -p $prefix/qc/fullres_trf

    touch ${prefix}/${moving_img.baseName}_registered.ome.tif
    touch ${prefix}/qc/refchns/elastix.log
    touch ${prefix}/qc/fullres_trf/transform.txt
    """
}
