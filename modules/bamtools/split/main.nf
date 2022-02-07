process BAMTOOLS_SPLIT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bamtools=2.5.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamtools:2.5.1--h9a82719_9' :
        'quay.io/biocontainers/bamtools:2.5.1--h9a82719_9' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bamtools \\
        split \\
        -in $bam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtools: \$( bamtools --version | grep -e 'bamtools' | sed 's/^.*bamtools //' )
    END_VERSIONS
    """
}
