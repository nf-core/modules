process CONCOCT_MERGECUTUPCLUSTERING {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/concoct:1.1.0--py312h245ed52_6':
        'biocontainers/concoct:1.1.0--py312h245ed52_6' }"

    input:
    tuple val(meta), path(clustering_csv)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$clustering_csv" == "${prefix}.csv") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    merge_cutup_clustering.py \\
        $args \\
        $clustering_csv \\
        > ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        concoct: \$(echo \$(concoct --version 2>&1) | sed 's/concoct //g' )
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$clustering_csv" == "${prefix}.csv") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        concoct: \$(echo \$(concoct --version 2>&1) | sed 's/concoct //g' )
    END_VERSIONS
    """
}
