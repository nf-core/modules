process BCFTOOLS_STATS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::bcftools=1.13' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.13--h3a49de5_0' :
        'quay.io/biocontainers/bcftools:1.13--h3a49de5_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*stats.txt"), emit: stats
    path  "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    bcftools stats $args $vcf > ${prefix}.bcftools_stats.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
