process HMMER_HMMPRESS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_1' :
        'quay.io/biocontainers/hmmer:3.4--hdbdd923_1' }"

    input:
    tuple val(meta), path(hmmfile)

    output:
    tuple val(meta), path("*.h3?"), emit: compressed_db
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), emit: versions_hmmer, topic: versions


    script:
    def args = task.ext.args ?: ''

    """
    hmmpress \\
        $args \\
        ${hmmfile}
    """

    stub:
    def prefix = task.ext.prefix ?: "stub"

    """
    touch ${prefix}.h3m
    touch ${prefix}.h3i
    touch ${prefix}.h3f
    touch ${prefix}.h3p
    """
}
