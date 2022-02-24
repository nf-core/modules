process BISCUIT_BSCONV {
    tag "$meta.id"
    label 'process_long'

    conda (params.enable_conda ? "bioconda::biscuit=1.0.2.20220113" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biscuit:1.0.2.20220113--h81a5ba2_0':
        'quay.io/biocontainers/biscuit:1.0.2.20220113--h81a5ba2_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(index)

    output:
    tuple val(meta), path("*.bsconv.bam"), emit: bsconv_bam
    path "versions.yml"                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    INDEX=`find -L ./ -name "*.bis.amb" | sed 's/.bis.amb//'`

    biscuit bsconv \\
        $args \\
        \$INDEX \\
        $bam \\
        ${prefix}.bsconv.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$(echo \$(biscuit version 2>&1) | sed 's/^.*BISCUIT Version: //; s/Using.*\$//')
    END_VERSIONS
    """
}
