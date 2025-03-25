process MELON {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/melon:0.2.5--pyhdfd78af_0':
        'biocontainers/melon:0.2.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path(database)

    output:
    tuple val(meta), path("**/*.tsv")    , emit: tsv_output
    tuple val(meta), path("**/*.json")   , emit: json_output
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    melon \\
        $reads \\
        --db $database \\
        --output ${prefix} \\
        --threads $task.cpus \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        melon: \$(melon -v)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.tsv
    touch ${prefix}/${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        melon: \$(melon -v)
    END_VERSIONS
    """
}
