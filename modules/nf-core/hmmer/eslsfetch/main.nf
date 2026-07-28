process HMMER_ESLSFETCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hb6cb901_4' :
        'quay.io/biocontainers/hmmer:3.4--hb6cb901_4' }"

    input:
    tuple val(meta), path(seqfile), path(keyfile), path(ssi)

    output:
    tuple val(meta), path("*.ssi", includeInputs: true), emit: ssi
    tuple val(meta), path("${prefix}.fasta"), emit: sequences, optional: true
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), emit: versions_hmmer, topic: versions
    tuple val("${task.process}"), val('easel'), eval("esl-sfetch -h | sed '2!d;s/^# Easel *//;s/ .*//'"), emit: versions_easel, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    def index_command = ssi ? '' : "esl-sfetch --index $seqfile"
    def fetch_command = keyfile ? "esl-sfetch -f ${args} -o ${prefix}.fasta ${seqfile} ${keyfile}" : ''
    """
    $index_command

    $fetch_command
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def index_command = ssi ? '' : "touch ${seqfile}.ssi"
    def fetch_command = keyfile ? "touch ${prefix}.fasta" : ''
    """
    $index_command
    $fetch_command
    """
}
