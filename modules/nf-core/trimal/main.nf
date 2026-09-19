process TRIMAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimal:1.5.0--h9948957_2':
        'quay.io/biocontainers/trimal:1.5.0--h9948957_2' }"

    input:
    tuple val(meta), path(aln)
    val out_format

    output:
    tuple val(meta), path("${prefix}.${out_extension}"), emit: trimal
    tuple val(meta), path("${prefix}.html")            , emit: summary, optional: true
    tuple val("${task.process}"), val('trimal'), eval("trimal --version | sed '2!d;s/trimAl v//;s/ .*//'"), emit: versions_trimal, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    out_extension = out_format ? out_format : "trimal"
    """
    trimal \\
        -in ${aln} \\
        -out ${prefix}.${out_extension} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    out_extension = out_format ? out_format : "trimal"
    """
    touch ${prefix}.${out_extension}
    touch ${prefix}.html
    """
}
