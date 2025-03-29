process BWA_ALN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bf/bf7890f8d4e38a7586581cb7fa13401b7af1582f21d94eef969df4cea852b6da/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:56c9f8d5201889a4' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(index)

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
        INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

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
        INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

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
