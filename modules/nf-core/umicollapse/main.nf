process UMICOLLAPSE {
    tag "$meta.id"
    label "process_medium"

    conda "bioconda::umicollapse=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umicollapse:1.0.0--hdfd78af_1' :
        'biocontainers/umicollapse:1.0.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: log
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    umicollapse \\
        bam \\
        -i $bam \\
        -o ${prefix}.bam \\
        $args

    mv .command.log ${prefix}_UMICollapse.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umicollapse: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dedup.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umicollapse: $VERSION
    END_VERSIONS
    """
}
