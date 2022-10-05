process BWA_ALN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bwa=0.7.17" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--h5bf99c6_8' :
        'quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8' }"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.sai"), emit: sai
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

        bwa aln \\
            $args \\
            -t $task.cpus \\
            -f ${prefix}.sai \\
            \$INDEX \\
            ${reads}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        END_VERSIONS
        """
    } else {
        """
        INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

        bwa aln \\
            $args \\
            -t $task.cpus \\
            -f ${prefix}.1.sai \\
            \$INDEX \\
            ${reads[0]}

        bwa aln \\
            $args \\
            -t $task.cpus \\
            -f ${prefix}.2.sai \\
            \$INDEX \\
            ${reads[1]}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        END_VERSIONS
        """
    }
}
