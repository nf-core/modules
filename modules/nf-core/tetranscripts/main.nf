process TETRANSCRIPTS {
    tag "$meta_c.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tetranscripts:2.2.3--pyh7cba7a3_0':
        'biocontainers/tetranscripts:2.2.3--pyh7cba7a3_0' }"

    input:
    tuple val(meta_t), path(bam_t)
    tuple val(meta_c), path(bam_c)
    tuple val(meta_ggtf), path(g_gtf)
    tuple val(meta_tegtf), path(te_gtf)

    output:
    tuple val(meta_t), path("*.cntTable"), emit: countTable
    tuple val(meta_t), path("*.R"), emit: log2fc
    tuple val(meta_t), path("*_analysis.txt"), emit: analysis, optional: true
    tuple val(meta_t), path("*_gene_TE.txt"), emit: sigdiff, optional: true
    tuple val("${task.process}"), val('tetranscripts'), eval("tetranscripts version | sed '1!d;s/.* //'"), emit: versions_tetranscripts, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta_c.id}"
// Join multiple BAM files with spaces for -t and -c arguments
    def treatment_bams = [bam_t].flatten().join(' ')
    def control_bams = [bam_c].flatten().join(' ')
    """
    TEtranscripts \\
	-t ${treatment_bams} \\
	-c ${control_bams} \\
	--GTF $g_gtf \\
	--TE $te_gtf \\
	--project ${prefix} \\
        $args

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta_c.id}"
    """
    echo $args

    touch ${prefix}.R
    touch ${prefix}.cntTable
    touch ${prefix}_gene_TE_analysis.txt
    touch ${prefix}_sigdiff_gene_TE.txt

    """
}
