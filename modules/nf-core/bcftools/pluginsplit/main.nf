process BCFTOOLS_PLUGINSPLIT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path(samples)
    path(groups)
    path(regions)
    path(targets)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}")  , emit: vcf
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def samples_arg = samples ? "--samples-file ${samples}" : ""
    def groups_arg  = groups  ? "--groups-file ${groups}"   : ""
    def regions_arg = regions ? "--regions-file ${regions}" : ""
    def targets_arg = targets ? "--targets-file ${targets}" : ""

    """
    bcftools plugin split \\
        ${vcf} \\
        ${samples_arg} \\
        ${groups_arg} \\
        ${regions_arg} \\
        ${targets_arg} \\
        --output ${prefix}

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

    def determination_file = samples ?: targets
    """
    cut -f 3 ${determination_file} | sed -e 's/\$/.${extension}/' > files.txt

    touch \$(<files.txt)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
