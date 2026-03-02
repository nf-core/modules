process WISECONDORX_NEWREF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wisecondorx:1.2.9--pyhdfd78af_0':
        'biocontainers/wisecondorx:1.2.9--pyhdfd78af_0' }"

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
