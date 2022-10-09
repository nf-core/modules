process ISMAPPER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ismapper=2.0.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ismapper:2.0.2--pyhdfd78af_1' :
        'quay.io/biocontainers/ismapper:2.0.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(reads), path(reference), path(query)

    output:
    tuple val(meta), path("results/*"), emit: results
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ismap \\
        $args \\
        --t $task.cpus \\
        --output_dir results \\
        --queries $query \\
        --reference $reference \\
        --reads $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ismapper: \$( echo \$( ismap --version 2>&1 ) | sed 's/^.*ismap //' )
    END_VERSIONS
    """
}
