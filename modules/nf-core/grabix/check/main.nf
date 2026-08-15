process GRABIX_CHECK {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grabix:0.1.8--hdcf5f25_9':
        'quay.io/biocontainers/grabix:0.1.8--hdcf5f25_9' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), stdout, emit: compress_bgzip
    tuple val("${task.process}"), val('grabix'), eval("grabix |& sed -n 's/^version: //p'"), topic: versions, emit: versions_grabix

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    grabix check ${input} | tr -d '\\n'
    """

    stub:

    """
    echo yes | tr -d '\\n'
    """
}
