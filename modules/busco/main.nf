process BUSCO {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::busco=5.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(fastas)
    path(lineage)

    output:
    tuple val(meta), path("short_summary.*.txt"),    emit: short_summary
    tuple val(meta), path("full_table.tsv"),         emit: full_table
    tuple val(meta), path("missing_busco_list.tsv"), emit: missing_busco_list
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    busco \\
        -i \\
        -l \\
        -o \\
        -c $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
