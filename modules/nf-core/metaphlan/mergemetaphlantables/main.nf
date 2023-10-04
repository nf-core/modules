process METAPHLAN_MERGEMETAPHLANTABLES {
    label 'process_single'

    conda "bioconda::metaphlan=4.0.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metaphlan:4.0.6--pyhca03a8a_0' :
        'biocontainers/metaphlan:4.0.6--pyhca03a8a_0' }"

    input:
    tuple val(meta), path(profiles)

    output:
    tuple val(meta), path("${prefix}.txt") , emit: txt
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_metaphlan_tables.py \\
        $args \\
        -o ${prefix}.txt \\
        ${profiles}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(metaphlan --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """
}
