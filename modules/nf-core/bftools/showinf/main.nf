process BFTOOLS_SHOWINF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bftools:8.0.0--hdfd78af_0':
        'biocontainers/bftools:8.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*.xml.gz")                                                             , emit: xml
    tuple val("${task.process}"), val("bftools"), eval("showinf -version | sed -n '1s/[^ ]* //p'"), emit: versions_bftools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    export BF_FLAGS='-XX:+PerfDisableSharedMem'
    showinf -nopix -no-upgrade -omexml-only \\
        $args \\
        $image | gzip > ${prefix}.xml.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo '<?xml version="1.0" encoding="UTF-8">' | gzip > ${prefix}.xml.gz
    """
}
