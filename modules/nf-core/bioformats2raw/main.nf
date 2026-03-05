process BIOFORMATS2RAW {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/523f11e0352a16c92c05b2baaccfc8e42c7b4d47e4997820478a90f45efa2fbc/data' :
        'community.wave.seqera.io/library/bioformats2raw:0.9.4--3eec45888b3759e5'}"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*ome.zarr"),         emit: omezarr
    tuple val("${task.process}"), val('bioformats2raw'), eval('bioformats2raw --version |& sed -n "1s/Version = //p"')         , emit: versions_bioformats2raw, topic: versions
    tuple val("${task.process}"), val('bio-formats'), eval('bioformats2raw --version |& sed -n "2s/Bio-Formats version = //p"'), emit: versions_bioformats, topic: versions
    tuple val("${task.process}"), val('ngff'), eval('bioformats2raw --version |& sed -n "3s/NGFF specification version = //p"'), emit: versions_ngff, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bioformats2raw \\
        $image \\
        ${prefix}.ome.zarr \\
        --max-workers $task.cpus \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir ${prefix}.ome.zarr
    """
}
