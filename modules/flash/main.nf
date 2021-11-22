process FLASH {
    tag "$meta.id"
    label 'process_medium'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/flash:1.2.11--hed695b0_5"
    } else {
        container "quay.io/biocontainers/flash:1.2.11--hed695b0_5"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml"                , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    flash \\
        $args \\
        -o ${prefix} \\
        -z \\
        ${reads[0]} \\
        ${reads[1]}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(flash --version 2>&1) | sed 's/^.*FLASH v//; s/ .*\$//')
    END_VERSIONS
    """
}
