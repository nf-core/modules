process SEACR_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::seacr=1.3 conda-forge::r-base=4.0.2 bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-03bfeb32fe80910c231f630d4262b83677c8c0f4:f4bb19b68e66de27e4c64306f951d5ff11919931-0' :
        'biocontainers/mulled-v2-03bfeb32fe80910c231f630d4262b83677c8c0f4:f4bb19b68e66de27e4c64306f951d5ff11919931-0' }"

    input:
    tuple val(meta), path(bedgraph), path(ctrlbedgraph)
    val (threshold)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def function_switch = ctrlbedgraph ? "$ctrlbedgraph" : "$threshold"
    def VERSION = '1.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    SEACR_1.3.sh \\
        $bedgraph \\
        $function_switch \\
        $args \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seacr: $VERSION
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
