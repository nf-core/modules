process BBMAP_INDEX {
    tag "$fasta"
    label 'process_long'

    conda (params.enable_conda ? "bioconda::bbmap=38.92" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:38.92--he522d1c_0' :
        'quay.io/biocontainers/bbmap:38.92--he522d1c_0' }"

    input:
    path fasta

    output:
    path 'ref'                    , emit: index
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bbmap.sh \\
        ref=${fasta} \\
        $args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh)
    END_VERSIONS
    """
}
