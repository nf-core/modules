// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SMOOTHXG {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4b/4b9fdcd26aa8b82cdf91b6d9afa9b092be7d5eabc810b462e410e6b501b53ce1/data'
        : 'community.wave.seqera.io/library/smoothxg:0.8.2--49c8fc69d4748df5'}"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*smoothxg.gfa"), emit: gfa
    path("*.maf") , optional: true, emit: maf
    tuple val("${task.process}"), val('smoothxg'), eval("smoothxg --version 2>&1 | sed 's/^v//; s/-.*//'"  ), topic: versions, emit: versions_smoothxg

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    smoothxg \\
        --threads=$task.cpus \\
        --gfa-in=${gfa} \\
        --smoothed-out=${prefix}.smoothxg.gfa \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.smoothxg.gfa
    """
}
