
process TRUVARI_CONSISTENCY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/truvari:5.4.0--pyhdfd78af_0':
        'biocontainers/truvari:5.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("*.{txt,json}") , emit: consistency
    tuple val("${task.process}"), val('truvari'), eval("truvari version | sed 's/Truvari v//'"), topic: versions, emit: versions_truvari

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("-j") ? "json" : "txt"

    """
    truvari \\
        consistency \\
        $args \\
        $vcfs > ${prefix}.${extension}

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("-j") ? "json" : "txt"

    """
    touch ${prefix}.${extension}
    """
}
