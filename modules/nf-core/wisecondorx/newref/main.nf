process WISECONDORX_NEWREF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/10a1cb134d692ea5d26cb7b4a5a2e83fe28eacac5b637a8eb6ca4d9532602222/data':
        'community.wave.seqera.io/library/wisecondorx:1.3.2--7ece9fb804446823' }"

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
