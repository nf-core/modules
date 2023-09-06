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
    tuple val(meta), path("${prefix}.bam") , emit: primary_filtered_bam
    tuple val(meta), path("${prefix2}.bam"), emit: contamination_bam
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_primary"
    prefix2 = task.ext.prefix2 ?: "${meta.id}_contaminant"

    if ("$primary_aligned_bam" == "${prefix}.bam"  | "$contaminant_aligned_bam" == "${prefix}.bam"  )
        error "Input and output names for the primary-genome bam file are the same, use \"task.ext.prefix\" to disambiguate!"
    if ("$primary_aligned_bam" == "${prefix2}.bam" | "$contaminant_aligned_bam" == "${prefix2}.bam" )
        error "Input and output names for the contaminant-genome bam file are the same, use \"task.ext.prefix2\" to disambiguate!"
    if ("primary_aligned_bam"    == "contaminant_aligned_bam" )
        error "Input file names for the two bam files are the same, ensure they are renamed upstream."
    if ("${prefix}.bam"    == "${prefix2}.bam" )
        error "Output names for the two bam files are identical, use \"task.ext.prefix\" and \"task.ext.prefix2\" to disambiguate!"

    def VERSION = '2.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    bamcmp \\
        -1 $primary_aligned_bam \\
        -2 $contaminant_aligned_bam \\
        -A ${prefix}.bam \\
        -B ${prefix2}.bam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamcmp: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_primary"
    prefix2 = task.ext.prefix2 ?: "${meta.id}_contaminant"

    if ("$primary_aligned_bam" == "${prefix}.bam"  | "$contaminant_aligned_bam" == "${prefix}.bam"  )
        error "Input and output names for the primary-genome bam file are the same, use \"task.ext.prefix\" to disambiguate!"
    if ("$primary_aligned_bam" == "${prefix2}.bam" | "$contaminant_aligned_bam" == "${prefix2}.bam" )
        error "Input and output names for the contaminant-genome bam file are the same, use \"task.ext.prefix2\" to disambiguate!"
    if ("primary_aligned_bam"    == "contaminant_aligned_bam" )
        error "Input file names for the two bam files are the same, ensure they are renamed upstream."
    if ("${prefix}.bam"    == "${prefix2}.bam" )
        error "Output names for the two bam files are identical, use \"task.ext.prefix\" and \"task.ext.prefix2\" to disambiguate!"
    """
    touch ${prefix}.bam
    touch ${prefix2}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g")
    END_VERSIONS
    """

}
