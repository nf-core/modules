process BAMUTIL_TRIMBAM {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bamutil=1.0.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamutil:1.0.15--h2e03b76_1' :
        'quay.io/biocontainers/bamutil:1.0.15--h2e03b76_1' }"

    input:
    tuple val(meta), path(bam), val(trim_left), val(trim_right)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bam \\
        trimBam \\
        $bam \\
        ${prefix}.bam \\
        $args \\
        -L $trim_left \\
        -R $trim_right

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamutil: \$( echo \$( bam trimBam 2>&1 ) | sed 's/^Version: //;s/;.*//' )
    END_VERSIONS
    """
}
