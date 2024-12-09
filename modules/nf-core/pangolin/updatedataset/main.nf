process PANGOLIN_UPDATEDATASET {
    tag "$dbname"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pangolin:4.3.1--pyhdfd78af_0':
        'biocontainers/pangolin:4.3.1--pyhdfd78af_0' }"

    input:
    val(dbname)

    output:
    tuple val(dbname), path("${prefix}"), emit: db
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${dbname}"
    """
    export XDG_CACHE_HOME=/tmp/.cache
    mkdir -p ${prefix}

    pangolin \\
        $args \\
        --update-data \\
        --threads $task.cpus \\
        --datadir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangolin: \$(pangolin --version | sed "s/pangolin //g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${dbname}"
    """
    mkdir -p ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangolin: \$(pangolin --version | sed "s/pangolin //g")
    END_VERSIONS
    """
}
