process PRETEXTGRAPH {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/28/28a9e6732fa1e04dd010a1596648055ec4b5bd300a16fd2f7b2c10e6594ad787/data'
        : 'community.wave.seqera.io/library/pretextgraph:0.0.9--c1c46eb0f0de5e85' }"

    input:
    tuple val(meta), path(pretext), path(bedgraph), val(graph_name)

    output:
    tuple val(meta), path("*.pretext"), emit: pretext
    tuple val("${task.process}"), val('PretextGraph'), eval('PretextGraph 2>&1 | sed "/Version/!d; s/.*Version //"'), emit: versions_pretextgraph, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${graph_name}"
    if ("${prefix}.pretext" == "${pretext}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    zcat -f "${bedgraph}" | PretextGraph \\
        ${args} \\
        -i "${pretext}" \\
        -n "${graph_name}" \\
        -o "${prefix}.pretext"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${graph_name}"
    if ("${prefix}.pretext" == "${pretext}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    touch "${prefix}.pretext"
    """
}
