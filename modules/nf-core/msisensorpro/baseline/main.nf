process MSISENSORPRO_BASELINE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msisensor-pro:1.2.0--hfc31af2_0' :
        'biocontainers/msisensor-pro:1.2.0--hfc31af2_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(msisensor_scan)
    tuple val(meta3), path(configure)

    output:
    tuple val(meta), path("*_baseline.list"), emit: baseline
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    msisensor-pro \\
        baseline \\
        -d ${msisensor_scan} \\
        -i ${configure}
        -o $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor-pro: \$(msisensor-pro 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p')
    END_VERSIONS
    """
}
