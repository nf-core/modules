process BACKSUB {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ab/ab42ea4815af3fd8a6ab667a7cd682b17808bac13214ddab2467dbcb80ed5dd8/data' :
        'community.wave.seqera.io/library/pip_backsub:1f25bcdea7241dd7'}"

    input:
    tuple val(meta) , path(image)
    tuple val(meta2), path(markerfile)

    output:
    tuple val(meta), path("*.ome.tif"), emit: backsub_tif
    tuple val(meta2), path("*.csv")   , emit: markerout
    tuple val("${task.process}"), val('backsub'), eval('backsub --version'), emit: versions_backsub, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$image" == "${prefix}.ome.tif") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    backsub \
        --input $image \
        --markers $markerfile \
        --output "${prefix}.ome.tif" \
        --marker-output "${prefix}.csv" \
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.ome.tif"
    touch "${prefix}.csv"
    """
}
