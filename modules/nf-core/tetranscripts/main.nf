process TETRANSCRIPTS {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tetranscripts:2.2.3--pyh7cba7a3_0':
        'biocontainers/tetranscripts:2.2.3--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(bam_t)
    tuple val(meta_c), path(bam_c)
    tuple val(meta_ggtf), path(g_gtf)
    tuple val(meta_tegtf), path(te_gtf)

    output:
    // TODO nf-core: Update the information obtained from bio.tools and make sure that it is correct
    tuple val(meta), path("*.{cntTable}"), emit: countTable
    tuple val(meta_c), path("*.{R}"), emit: log2fc
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tetranscripts \\
	-t $bam_t \\
	-c $bam_c \\
	--GTF $g_gtf \\
	--TE $te_gtf \\
	--project ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tetranscripts: \$(tetranscripts --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}.R
    touch ${prefix}.cntTable

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tetranscripts: \$(tetranscripts --version)
    END_VERSIONS
    """
}
