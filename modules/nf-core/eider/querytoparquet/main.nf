process EIDER_QUERYTOPARQUET {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eider:0.1--hdfd78af_0' :
        'biocontainers/eider:0.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(sql)

    output:
    tuple val(meta), path("*.parquet"), emit: parquet
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    eider \
        $args \
        --skip-history \
        --verbose \
        --parameters prefix=$prefix \
        --query-path $sql

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eider: \$(eider --version 2>&1 | grep -o 'eider .*' | cut -f2 -d ' ')
    END_VERSIONS
    """
}
