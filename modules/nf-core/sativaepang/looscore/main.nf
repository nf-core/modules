process SATIVAEPANG_LOOSCORE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sativa-epang:0.9.3.4--py314hab16a5f_0' :
        'quay.io/biocontainers/sativa-epang:0.9.3.4--py314hab16a5f_0' }"

    input:
    tuple val(meta), path(refjson), path(taskdir)

    output:
    tuple val(meta), path("*.mis"), emit: mis
    tuple val("${task.process}"), val('sativaepang'), eval("grep -m1 -oE '[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+' \$(command -v sativa-epang)"), topic: versions, emit: versions_sativaepang

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Unlike looplace, loo-score only reads taskdir -- writes its own output via -o.
    """
    sativa-epang \\
        -r ${refjson} \\
        -n ${prefix} \\
        -o . \\
        -taskdir ${taskdir} \\
        -stage loo-score \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mis
    """
}
