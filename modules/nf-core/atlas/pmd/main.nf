process ATLAS_PMD {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::atlas=0.9.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/atlas:0.9.9--h082e891_0':
        'biocontainers/atlas:0.9.9--h082e891_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(pool_rg_txt)
    path(fasta)
    path(fai)

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
    def pool_rg_txt = pool_rg_txt ? "poolReadGroups=${pool_rg_txt}" : ""
    """
    atlas \\
        $pool_rg_txt \\
        task=PMD \\
        bam=${bam} \\
        fasta=${fasta} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        atlas: \$((atlas 2>&1) | grep Atlas | head -n 1 | sed -e 's/^[ \t]*Atlas //')
    END_VERSIONS
    """
}
