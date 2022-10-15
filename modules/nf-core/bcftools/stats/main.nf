process BCFTOOLS_STATS {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bcftools=1.15.1" : null)
    def container_image = "bcftools:1.15.1--h0ea216a_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(vcf), path(tbi)
    path regions
    path targets
    path samples

    output:
    tuple val(meta), path("*stats.txt"), emit: stats
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions_file = regions ? "--regions-file ${regions}" : ""
    def targets_file = targets ? "--targets-file ${targets}" : ""
    def samples_file =  samples ? "--samples-file ${samples}" : ""
    """
    bcftools stats \\
        $args \\
        $regions_file \\
        $targets_file \\
        $samples_file \\
        $vcf > ${prefix}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
