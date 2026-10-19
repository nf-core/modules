process TIDK_PLOT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9fc2352a02fcae6bc8ae7ba109692168150451bf188a98c1032545f85a0b21e2/data'
        : 'community.wave.seqera.io/library/tidk:0.2.7--e781b1e8244cc913'}"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.svg"), emit: svg
    tuple val("${task.process}"), val('tidk'), eval("tidk --version | sed 's/tidk //'"), emit: versions_tidk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tidk \\
        plot \\
        --output $prefix \\
        $args \\
        --tsv "$tsv"

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.svg

    """
}
