process KRAKEN2_BUILD {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f8c4015c836dd3ce5c118cfed97ec8259bab9e9d:2e0b144854b4a3d69b5df7a0340a60db846cc8bf-0':
        'biocontainers/mulled-v2-f8c4015c836dd3ce5c118cfed97ec8259bab9e9d:2e0b144854b4a3d69b5df7a0340a60db846cc8bf-0' }"

    input:
    tuple val(meta), path(db)
    val cleaning

    output:
    tuple val(meta), path("$prefix"), emit: db
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    runclean = cleaning ? "kraken2-build --clean --db ${db}" : ""
    """
    kraken2-build \\
        --build \\
        $args \\
        --threads ${task.cpus} \\
        --db ${db}
    $runclean 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "$prefix"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}
