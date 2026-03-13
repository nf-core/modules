process GCTA_MAKEBKSPARSE {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gcta:1.94.1--h9ee0642_0' :
        'biocontainers/gcta:1.94.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(grm_id), path(grm_bin), path(grm_n_bin)
    val(cutoff)

    output:
    tuple val(meta), path("${meta.id}_sp.grm.id"), path("${meta.id}_sp.grm.sp"), emit: sparse_grm_files
    tuple val("${task.process}"), val("gcta"), eval("gcta --version 2>&1 | head -n 1"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    gcta \\
        --grm ${meta.id} \\
        --make-bK-sparse ${cutoff} \\
        --out ${meta.id}_sp \\
        --thread-num ${task.cpus} \\
        ${args}
    """

    stub:
    """
    touch ${meta.id}_sp.grm.id
    touch ${meta.id}_sp.grm.sp
    """
}
