process MYKROBE_PREDICT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mykrobe:0.11.0--py39h2add14b_1':
        'biocontainers/mykrobe:0.11.0--py39h2add14b_1' }"

    input:
    tuple val(meta), path(seqs)
    val species

    output:
    tuple val(meta), path("${prefix}.csv") , emit: csv
    tuple val(meta), path("${prefix}.json"), emit: json
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mykrobe \\
        predict \\
        $args \\
        --species $species \\
        --threads $task.cpus \\
        --sample $prefix \\
        --format json_and_csv \\
        --output ${prefix} \\
        --seq $seqs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mykrobe: \$(echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mykrobe: \$(echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v//' )
    END_VERSIONS
    """
}
