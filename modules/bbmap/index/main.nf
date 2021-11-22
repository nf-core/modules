process BBMAP_INDEX {
    tag "$fasta"
    label 'process_long'

    conda (params.enable_conda ? "bioconda::bbmap=38.92" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bbmap:38.92--he522d1c_0"
    } else {
        container "quay.io/biocontainers/bbmap:38.92--he522d1c_0"
    }

    input:
    path fasta

    output:
    path 'ref'                    , emit: index
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    bbmap.sh \\
        ref=${fasta} \\
        $args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(bbversion.sh)
    END_VERSIONS
    """
}
