process HMMER_ESLALIPID {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hb6cb901_4' :
        'quay.io/biocontainers/hmmer:3.4--hb6cb901_4' }"

    input:
    tuple val(meta), path(aln)

    output:
    tuple val(meta), path("${prefix}.alipid.txt"), emit: pid
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), topic: versions, emit: versions_hmmer
    tuple val("${task.process}"), val('easel'), eval("esl-alipid -h | sed '2!d;s/^# Easel *//;s/ .*//'"), topic: versions, emit: versions_easel

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    esl-alipid \\
        $args \\
        $aln \\
        > ${prefix}.alipid.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alipid.txt
    """
}
