process WISECONDORX_NEWREF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/13/13af39819608398807612090d4b8af7dedb8db403967e71af22dbbeeb502ead1/data':
        'community.wave.seqera.io/library/wisecondorx:1.3.0--835c946afbce9082' }"

    input:
    tuple val(meta), path(inputs)

    output:
    tuple val(meta), path("*.npz"), emit: npz
    tuple val("${task.process}"), val('wisecondorx'), eval("pip list |& sed -n 's/wisecondorx *//p'"), emit: versions_wisecondorx, topic: versions

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
