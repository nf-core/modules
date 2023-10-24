process BCFTOOLS_PLUGINSCATTER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val(sites_per_chunk)
    val(scatter)
    path(scatter_file)
    path(regions)
    path(targets)

    output:
    tuple val(meta), path("*{vcf,vcf.gz,bcf,bcf.gz}")   , emit: scatter
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def mandatory_arg = sites_per_chunk ? "--nsites-per-chunk ${sites_per_chunk}" : scatter ? "--scatter ${scatter}" : "--scatter-file ${scatter_file}"
    def regions_arg = regions ? "--regions-file ${regions}" : ""
    def targets_arg = targets ? "--targets-file ${targets}" : ""
    """
    bcftools plugin scatter \\
        ${vcf} \\
        ${mandatory_arg} \\
        ${regions_arg} \\
        ${targets_arg} \\
        --output ${prefix} \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    mv ${prefix}/* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                "vcf"

    """
    touch ${prefix}1.${extension}
    touch ${prefix}2.${extension}
    touch ${prefix}3.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
