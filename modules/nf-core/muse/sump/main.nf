process MUSE_SUMP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(muse_call_txt)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '' // hands -G for WGS data and -E for WES data
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seqtype =
    def VERSION = '2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    MuSE \\
        sump \\
        $args \\
        -I $muse_call_txt \\
        -O ${prefix}.vcf \\
        -n $task.cpus    \\
        -D $reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        muse: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        muse: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
