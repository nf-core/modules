process PICARD_SCATTERINTERVALSBYNS {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.2.0--hdfd78af_0' :
        'biocontainers/picard:3.2.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta3), path(dict)

    output:
    tuple val(meta), path("*.interval_list"), emit: intervals
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    picard ScatterIntervalsByNs \\
        --REFERENCE ${fasta} \\
        --OUTPUT ${prefix}.interval_list \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard ScatterIntervalsByNs --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.interval_list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard ScatterIntervalsByNs --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
