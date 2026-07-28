process ICOUNTMINI_PEAKS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/icount-mini:3.0.1--pyh7cba7a3_0':
        'quay.io/biocontainers/icount-mini:3.0.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(bed), path(sigxls)

    output:
    tuple val(meta), path("*.peaks.bed.gz"), emit: peaks
    tuple val("${task.process}"), val('iCount-Mini'), eval("iCount-Mini -v"), emit: versions_icount_mini, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    iCount-Mini peaks \\
        ${bed} \\
        ${sigxls} \\
        ${prefix}.peaks.bed.gz \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.peaks.bed.gz
    """
}
