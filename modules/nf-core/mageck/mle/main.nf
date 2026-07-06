process MAGECK_MLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mageck:0.5.9.5--py39h1f90b4d_3':
        'quay.io/biocontainers/mageck:0.5.9.5--py39h1f90b4d_3' }"

    input:
    tuple val(meta), path(count_table)
    path(design_matrix)

    output:
    tuple val(meta), path("*.gene_summary.txt") , emit: gene_summary
    tuple val(meta), path("*.sgrna_summary.txt"), emit: sgrna_summary
    tuple val("${task.process}"), val("mageck"), eval("mageck -v"), emit: versions_mageck, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mageck \\
        mle \\
        ${args} \\
        --threads ${task.cpus} \\
        -k ${count_table} \\
        -d ${design_matrix} \\
        -n ${prefix}
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.gene_summary.txt
    touch ${prefix}.sgrna_summary.txt
    """


}
