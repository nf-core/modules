process HMMER_ESLSFETCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hb6cb901_4' :
        'quay.io/biocontainers/hmmer:3.4--hb6cb901_4' }"

    input:
    // seqfile must NOT be gzip-compressed: esl-sfetch --index cannot index a .gz file.
    // keyfile is optional (pass [] to skip): one name/accession per line for -f retrieval,
    // or a subsequence coordinate file for -Cf retrieval (set ext.args = '-C' to switch mode).
    tuple val(meta), path(seqfile), path(keyfile)

    output:
    tuple val(meta), path("*.ssi")          , emit: ssi
    tuple val(meta), path("${prefix}.fasta"), emit: sequences, optional: true
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), emit: versions_hmmer, topic: versions
    tuple val("${task.process}"), val('easel'), eval("esl-sfetch -h | sed '2!d;s/^# Easel *//;s/ .*//'"), emit: versions_easel, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    def fetch_command = keyfile ? "esl-sfetch -f ${args} -o ${prefix}.fasta ${seqfile} ${keyfile}" : ''
    """
    esl-sfetch --index $seqfile

    $fetch_command
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def fetch_command = keyfile ? "touch ${prefix}.fasta" : ''
    """
    touch ${seqfile}.ssi
    $fetch_command
    """
}
