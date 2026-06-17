process ICOUNTMINI_SIGXLS {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/icount-mini:3.0.1--pyh7cba7a3_0':
        'quay.io/biocontainers/icount-mini:3.0.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(bed)
    path segmentation

    output:
    tuple val(meta), path("*.sigxls.bed.gz"), emit: sigxls
    tuple val(meta), path("*.scores.tsv")   , emit: scores
    tuple val("${task.process}"), val('iCount-Mini'), eval("iCount-Mini -v"), emit: versions_icount_mini, topic: versions


    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    iCount-Mini sigxls \\
        ${segmentation} \\
        ${bed} \\
        ${prefix}.sigxls.bed.gz \\
        --scores ${prefix}.scores.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.sigxls.bed.gz
    touch ${prefix}.scores.tsv
    """
}
