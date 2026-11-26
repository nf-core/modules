// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process WISECONDORX_NEWREF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/41/4103a257615e467895ebb0ff8aca0c17f80e852466feb36736236a4ca2597abf/data'
        : 'community.wave.seqera.io/library/wisecondorx:1.3.3--dd93194eb6a11e3b'}"

    input:
    tuple val(meta), path(inputs)

    output:
    tuple val(meta), path("*.npz"), emit: npz
    tuple val("${task.process}"), val('wisecondorx'), eval("python -c \"import wisecondorx; print(wisecondorx.__version__)\""), emit: versions_wisecondorx, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    inputs.each { input ->
        if("${input}" == "${prefix}.npz") error "${input} has the same name as the output file, set prefix in module configuration to disambiguate!"
    }

    """
    WisecondorX \\
        newref \\
        *.npz \\
        ${prefix}.npz \\
        --cpus ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    inputs.each { input ->
        if("${input}" == "${prefix}.npz") error "${input} has the same name as the output file, set prefix in module configuration to disambiguate!"
    }

    """
    touch ${prefix}.npz
    """
}
