process BBMAP_SENDSKETCH {
    tag "$fasta"
    label 'process_high'

    conda ("bioconda::bbmap=39.01")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
        'quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    path  fasta

    output:
    path 'results'                , emit: results 
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """    
    sendsketch.sh \\
        in=${fasta} \\
        $args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
