process DRAGMAP_HASHTABLE {
    tag "$fasta"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::dragmap=1.2.1" : null)
        'https://depot.galaxyproject.org/singularity/dragmap:1.2.1--hd4ca14e_0':
        "${params.docker_registry ?: 'quay.io/biocontainers'}/dragmap:1.2.1--hd4ca14e_0" }

    input:
    path fasta

    output:
    path "dragmap"      , emit: hashmap
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir dragmap
    dragen-os \\
        --build-hash-table true \\
        --ht-reference $fasta \\
        --output-directory dragmap \\
        $args \\
        --ht-num-threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragmap: \$(echo \$(dragen-os --version 2>&1))
    END_VERSIONS
    """
}
