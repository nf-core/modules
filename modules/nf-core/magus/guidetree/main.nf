process MAGUS_GUIDETREE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magus-msa:0.2.0--pyhdfd78af_0':
        'quay.io/biocontainers/magus-msa:0.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tree"), emit: tree
    tuple val("${task.process}"), val('magus'), eval('magus --version'), emit: versions_magus, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    magus \\
        -np ${task.cpus} \\
        -i ${fasta} \\
        -o ${prefix}.tree \\
        --onlyguidetree TRUE \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tree
    """
}
