process CONIFER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::conifer=1.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/conifer%3A1.0.2--he4a0461_0':
        'biocontainers/conifer:1.0.2--he4a0461_0' }"

    input:
    tuple val(meta), path(kraken_result)
    path kraken_taxon_db
    
    output:
    tuple val(meta), path("*.score"), emit: score
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    conifer \\
        $args \\
        --input $kraken_result \\
        --db $kraken_taxon_db > ${prefix}.score

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        conifer: \$(echo \$(conifer --version 2>&1) | sed 's/^.*Conifer //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.score

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        conifer: \$(echo \$(conifer --version 2>&1) | sed 's/^.*Conifer //')
    END_VERSIONS
    """
}