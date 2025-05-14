process FASTQUTILS_INFO {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq_utils:0.25.2--h96c455f_2':
        'biocontainers/fastq_utils:0.25.2--h96c455f_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.txt"), emit: info
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.25.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    fastq_info \\
        $args \\
        ${reads}

    echo "fastq_utils fastq_info ran and found no issues with ${reads}" > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq_utils: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.25.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqutils: ${VERSION}
    END_VERSIONS
    """
}
