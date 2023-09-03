process BAMCMP {
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::bamcmp=2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamcmp:2.2--h05f6578_0' :
        'biocontainers/bamcmp:2.2--h05f6578_0' }"

    input:
    tuple val(meta), path(primary_aligned_bam), path(contaminant_aligned_bam)

    output:
    tuple val(meta), path("${$prefix1}.bam"), emit: primary_filtered_bam
    tuple val(meta), path("${$prefix2}.bam"), emit: contamination_bam
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def prefix1 = task.ext.prefix1 ?: "${meta.prefix}_primary"
    def prefix2 = task.ext.prefix2 ?: "${meta.prefix}_contaminant"
    if ("$primary_aligned_bam" == "${prefix1}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    if ("$contaminant_aligned_bam" == "${prefix1}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    if ("$primary_aligned_bam" == "${prefix2}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    if ("$contaminant_aligned_bam" == "${prefix2}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    def VERSION = '2.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    bamcmp \\
        -1 $primary_aligned_bam \\
        -2 $contaminant_aligned_bam \\
        -A ${prefix1}.bam \\
        -B ${prefix2}.bam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamcmp: $VERSION
    END_VERSIONS
    """

}
