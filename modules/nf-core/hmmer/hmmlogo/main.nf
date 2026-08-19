process HMMER_HMMLOGO {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hb6cb901_4' :
        'quay.io/biocontainers/hmmer:3.4--hb6cb901_4' }"

    input:
    tuple val(meta), path(hmm)

    output:
    tuple val(meta), path("${prefix}.hmmlogo.txt"),                          emit: logo
    tuple val("${task.process}"), val('hmmer'), eval("hmmlogo -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), topic: versions, emit: versions_hmmer

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    hmmlogo \\
        $args \\
        $hmm \\
        > ${prefix}.hmmlogo.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hmmlogo.txt
    """
}
