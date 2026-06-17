process HMMER_ESLREFORMAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hb6cb901_4' :
        'quay.io/biocontainers/hmmer:3.4--hb6cb901_4' }"

    input:
    tuple val(meta), path(seqfile)
    val postprocessing_script

    output:
    tuple val(meta), path("*.*.gz"), emit: seqreformated
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), emit: versions_hmmer, topic: versions
    tuple val("${task.process}"), val('easel'), eval("esl-reformat -h | sed '2!d;s/^# Easel *//;s/ .*//'"), emit: versions_easel, topic: versions


    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def suffix   = args ? args.trim().tokenize(" ")[-1] : "sequences"
    """
    esl-reformat \\
        $args \\
        $seqfile \\
        $postprocessing_script \\
        | gzip -c > ${prefix}.${suffix}.gz
    """

    stub:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def suffix   = args ? args.trim().tokenize(" ")[-1] : "sequences"

    """
    echo "" | gzip > ${prefix}.${suffix}.gz
    """
}
