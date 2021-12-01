process DRAGMAP_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::dragmap=1.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dragmap:1.2.1--hd4ca14e_0':
        'quay.io/biocontainers/dragmap:1.2.1--hd4ca14e_0' }"

    input:
    tuple val(meta), path(reads)
    path  hashmap

    output:
    tuple val(meta), path("*.sam"), emit: sam
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    if (meta.single_end) {
        """
        dragen-os \\
            -r $hashmap \\
            -1 $reads \\
            --output-directory . \\
            --output-file-prefix $prefix \\
            --num-threads $task.cpus \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        dragmap: \$(echo \$(dragen-os --version 2>&1))
        END_VERSIONS
        """
    } else {
        """
        dragen-os \\
            -r $hashmap \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            --output-directory . \\
            --output-file-prefix $prefix \\
            --num-threads $task.cpus \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        dragmap: \$(echo \$(dragen-os --version 2>&1))
        END_VERSIONS
        """
    }
}
