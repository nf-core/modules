process RP3NET_RP3 {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://rp3net/rp3net:v0.0.1' :
        'rp3net/rp3net:v0.0.1' }"

    input:
    tuple val(meta), path(fasta), path(checkpoint)

    output:
    tuple val(meta), path("${prefix}.csv"), emit: scores
    tuple val("${task.process}"), val('rp3net'), eval("pip show RP3Net 2>&1 | awk '/^Version/{print \$2}'"), topic: versions, emit: versions_rp3net

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    rp3 \\
        -p ${checkpoint} \\
        -f ${fasta} \\
        -o ${prefix}.csv \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    """
}
