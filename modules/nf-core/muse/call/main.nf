process MUSE_CALL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/muse_sump:6020175d1ed543c4':
        'community.wave.seqera.io/library/muse_sump:3847abd544ae3eb6' }"

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.MuSE.txt"), emit: txt
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    MuSE \\
        call \\
        $args \\
        -f $reference \\
        -O ${prefix}  \\
        -n $task.cpus \\
        $tumor_bam    \\
        $normal_bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MuSE: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.MuSE.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MuSE: ${VERSION}
    END_VERSIONS
    """
}
