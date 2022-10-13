process MYKROBE_PREDICT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mykrobe=0.11.0" : null)
    def container_image = "/mykrobe:0.11.0--py39h2add14b_1"
                                                  container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

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
}
