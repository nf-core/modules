process SLIMFASTQ {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::slimfastq=2.04"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/slimfastq:2.04--h87f3376_2':
        'biocontainers/slimfastq:2.04--h87f3376_2' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.sfq"), emit: sfq
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.04' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    if (meta.single_end) {
        """
        gzip -d -c '${fastq}' | slimfastq \\
            $args \\
            -f '${prefix}.sfq'

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            slimfastq: ${VERSION}
        END_VERSIONS
        """
    } else {
        """
        gzip -d -c '${fastq[0]}' | slimfastq \\
            $args \\
            -f '${prefix}_1.sfq'

        gzip -d -c '${fastq[1]}' | slimfastq \\
            $args \\
            -f '${prefix}_2.sfq'

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            slimfastq: ${VERSION}
        END_VERSIONS
        """
    }
}
