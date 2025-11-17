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
    tuple val(meta), path("*.xml"), emit: xml
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    export BF_FLAGS='-XX:+PerfDisableSharedMem'
    showinf -nopix -no-upgrade -omexml-only \\
        $args \\
        $image > ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        showinf: \$(showinf -version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo '<?xml version="1.0" encoding="UTF-8">' > ${prefix}.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        showinf: \$(showinf -version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """
}
