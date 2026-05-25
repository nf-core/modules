process GAPPA_EXAMINEGRAFT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gappa:0.8.0--h9a82719_0':
        'quay.io/biocontainers/gappa:0.8.0--h9a82719_0' }"

    input:
    tuple val(meta), path(jplace)

    output:
    tuple val(meta), path("*.newick"), emit: newick
    tuple val("${task.process}"), val('gappa'), eval("gappa --version 2>&1 | sed 's/v//'"), emit: versions_gappa, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gappa \\
        examine \\
        graft \\
        ${args} \\
        --threads ${task.cpus} \\
        --file-prefix ${prefix}. \\
        --jplace-path ${jplace}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.newick
    """
}
