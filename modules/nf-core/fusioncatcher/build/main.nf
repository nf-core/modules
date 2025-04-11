process FUSIONCATCHER_BUILD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fusioncatcher:1.33--hdfd78af_5':
        'biocontainers/fusioncatcher:1.33--hdfd78af_5' }"

    input:
    val(meta)

    output:
    tuple val(meta), path("${prefix}")  , emit: reference
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    fusioncatcher-build \\
        ${args} \\
        --output=${prefix} \\
        --threads=${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: "\$(fusioncatcher --version 2>&1 | awk '{print \$2}')"
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: "\$(fusioncatcher --version 2>&1 | awk '{print \$2}')"
    END_VERSIONS
    """
}
