process IVAR_TRIM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ivar:1.4.4--h077b44d_0' :
        'quay.io/biocontainers/ivar:1.4.4--h077b44d_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path('*.log'), emit: log
    tuple val("${task.process}"), val('ivar'), eval("ivar version | sed -n 's|iVar version \\(.*\\)|\\1|p'"), emit: versions_ivar, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ivar trim \\
        ${args} \\
        -i ${bam} \\
        -b ${bed} \\
        -p ${prefix} \\
        > ${prefix}.ivar.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ivar.log
    touch ${prefix}.bam
    """
}
