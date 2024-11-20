process REGTOOLS_JUNCTIONSEXTRACT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::regtools=0.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/regtools:0.5.0--he941832_0' :
        'biocontainers/regtools:0.5.0--he941832_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.junc"), emit: junc
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    regtools junctions extract \\
        $args \\
        -o ${prefix}.junc \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regtools: \$(regtools --version 2>&1 | sed 's/^.*regtools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.junc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
