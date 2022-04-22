process ATLAS_PMD {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::atlas=0.9.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai), path(read_group_settings)

    output:
    tuple val(meta), path("*_PMD_input_Empiric.txt")    , emit: empiric
    tuple val(meta), path("*_PMD_input_Exponential.txt"), emit: exponential
    tuple val(meta), path("*_PMD_Table_counts.txt")     , emit: counts
    tuple val(meta), path("*_PMD_Table.txt")            , emit: table
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_group_settings = read_group_settings ? "${read_group_settings}" : ""
    """
    atlas \\
        $read_group_settings \\
        task=PMD \\
        bam=$bam \\
        fasta=$fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        atlas: \$(echo \$(atlas 2>&1) | grep Atlas | head -n 1 | sed -e 's/^[ \t]*Atlas //')
    END_VERSIONS
    """
}
