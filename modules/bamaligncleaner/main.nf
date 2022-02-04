process BAMALIGNCLEANER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bamaligncleaner=0.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamaligncleaner:0.2.1--pyhdfd78af_0' :
        'quay.io/biocontainers/bamaligncleaner:0.2.1--pyhdfd78af_0' }"

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
    bamAlignCleaner \\
        $args \\
        -o ${prefix}.bam \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamaligncleaner: \$(bamAlignCleaner --version | sed 's/.*version //')
    END_VERSIONS
    """
}
