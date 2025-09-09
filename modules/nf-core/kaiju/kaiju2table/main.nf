process KAIJU_KAIJU2TABLE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kaiju:1.10.0--h43eeafb_0':
        'biocontainers/kaiju:1.10.0--h43eeafb_0' }"

    input:
    tuple val(meta), path(input)
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
    dbnames=`find -L ${db} -name "*names.dmp"`
    kaiju2table   $args \\
        -t \$dbnodes \\
        -n \$dbnames \\
        -r ${taxon_rank} \\
        -o ${prefix}.txt \\
        ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaiju: \$(echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaiju: \$(echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' ))
    END_VERSIONS
    """
}
