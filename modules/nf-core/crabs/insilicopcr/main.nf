process CRABS_INSILICOPCR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crabs:1.0.7--pyhdfd78af_0':
        'biocontainers/crabs:1.0.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(crabsdb)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: txt
    tuple val("${task.process}"), val('crabs'), eval("crabs --help 2>/dev/null | grep -oE 'v[0-9.]+' | cut -c2-"), emit: versions_crabs, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}.insilicopcr"
    """
    crabs --in-silico-pcr \\
        --input ${crabsdb} \\
        --output ${prefix}.txt \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.insilicopcr"
    """
    touch ${prefix}.txt
    """
}
