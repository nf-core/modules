process BIOFORMATS2RAW {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7d1466c4737fc34f9f3b63117bffda79092a7a59c795b14f4d49bdc9ebf8512d/data' :
        'community.wave.seqera.io/library/bioformats2raw:0.12.1--503439f3c2940fe1'}"

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
