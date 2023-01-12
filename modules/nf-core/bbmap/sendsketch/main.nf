process BBMAP_SENDSKETCH {
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bbmap=39.01" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.txt"), emit: results 
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """    
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    sendsketch.sh \\
        in=${fasta}.fasta \\
        outsketch=${prefix}.txt \\
        $args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
