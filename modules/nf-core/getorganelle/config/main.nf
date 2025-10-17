process GETORGANELLE_CONFIG {
    tag "${organelle_type}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/getorganelle:1.7.7.0--pyh7cba7a3_0':
        'biocontainers/getorganelle:1.7.7.0--pyh7cba7a3_0' }"

    input:
    val(organelle_type)

    output:
    tuple val(organelle_type), path("getorganelle"), emit: db
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    get_organelle_config.py \\
        $args \\
        -a ${organelle_type} \\
        --config-dir getorganelle

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorganelle: \$(get_organelle_config.py --version | sed 's/^GetOrganelle v//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    mkdir -p getorganelle/{LabelDatabase,SeedDatabase}
    touch getorganelle/{LabelDatabase,SeedDatabase}/${organelle_type}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorganelle: \$(get_organelle_config.py --version | sed 's/^GetOrganelle v//g')
    END_VERSIONS
    """
}
