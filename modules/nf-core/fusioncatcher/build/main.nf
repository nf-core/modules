process FUSIONCATCHER_BUILD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fusioncatcher:1.33--hdfd78af_5':
        'quay.io/biocontainers/fusioncatcher:1.33--hdfd78af_5' }"

    input:
    val(meta)

    output:
    tuple val(meta), path("${prefix}"), emit: reference
    tuple val("${task.process}"), val('fusioncatcher'), eval("fusioncatcher --version |& sed 's/.* //'"), topic: versions, emit: versions_fusioncatcher

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    fusioncatcher-build \\
        ${args} \\
        --output=${prefix} \\
        --threads=${task.cpus}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    """
}
