process KAIJU_KAIJU2TABLE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::kaiju=1.8.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kaiju:1.8.2--h5b5514e_1':
        'biocontainers/kaiju:1.8.2--h2e03b76_0' }"

    input:
    tuple val(meta), path(results)
    path db
    val taxon_rank

    output:
    tuple val(meta), path('*.txt'), emit: summary
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dbnodes=`find -L ${db} -name "*nodes.dmp"`
    dbname=`find -L ${db} -name "*.fmi" -not -name "._*"`
    kaiju2table   $args \\
        -t \$dbnodes \\
        -n \$dbname \\
        -r ${taxon_rank} \\
        -o ${prefix}.txt \\
        ${results}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaiju: \$(echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' ))
    END_VERSIONS
    """
}
