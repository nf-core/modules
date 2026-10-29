// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process HYPHY_GARD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed3d742a4b3281f09189e5010b1b26855e2bdc5259dc8bdf1121a7b1939e8ca6/data'
:         'community.wave.seqera.io/library/hyphy:2.5.101--da14606d9f2ca413' }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*_output.txt"), emit: log
    tuple val("${task.process}"), val('hyphy'), eval("hyphy --version | sed -E 's/.*HYPHY ([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/'"), topic: versions, emit: versions_hyphy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hyphy gard \\
        CPU=${task.cpus} \\
        --alignment ${alignment} \\
        --output ${prefix}_GARD.json \\
        ${args} \\
        > ${prefix}_GARD_output.txt
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}_GARD.json
    touch ${prefix}_GARD_output.txt
    """
}
