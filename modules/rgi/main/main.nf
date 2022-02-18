process RGI_MAIN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::rgi=5.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:5.2.1--pyha8f3691_2':
        'quay.io/biocontainers/rgi:5.2.1--pyha8f3691_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.txt") , emit: tsv
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rgi \\
        main \\
        $args \\
        --num_threads $task.cpus \\
        --output_file $prefix \\
        --input_sequence $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version)
    END_VERSIONS
    """
}
