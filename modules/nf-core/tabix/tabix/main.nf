process TABIX_TABIX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92859404d861ae01afb87e2b789aebc71c0ab546397af890c7df74e4ee22c8dd/data' :
        'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa' }"

    input:
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path("*.{tbi,csi}"), emit: index
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    tabix \\
        --threads $task.cpus \\
        $args \\
        $tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(tabix --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def index = args.contains("-C ") || args.contains("--csi") ? "csi" : "tbi"
    """
    touch ${tab}.${index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(tabix --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
