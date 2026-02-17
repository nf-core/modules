process MACSYFINDER_DOWNLOAD {
    tag "${model_name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macsyfinder:2.1.6--pyhdfd78af_0' :
        'biocontainers/macsyfinder:2.1.6--pyhdfd78af_0' }"

    input:
    val model_name

    output:
    path "models/*"     , emit: models
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # macsydata installs models into the current directory by default
    # We'll create a models directory to store them
    mkdir -p models

    macsydata install \\
        --target models \\
        ${args} \\
        ${model_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macsyfinder: \$(macsyfinder --version 2>&1 | sed 's/^.*MacSyFinder //; s/ .*\$//')
        macsydata: \$(macsydata --version 2>&1 | sed 's/^.*macsydata //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p models
    touch models/${model_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macsyfinder: \$(macsyfinder --version 2>&1 | sed 's/^.*MacSyFinder //; s/ .*\$//')
        macsydata: \$(macsydata --version 2>&1 | sed 's/^.*macsydata //; s/ .*\$//')
    END_VERSIONS
    """
}
