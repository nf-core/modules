nextflow.preview.types = true

process tetranscripts {
    tag "$meta_c.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tetranscripts:2.2.3--pyh7cba7a3_0':
        'biocontainers/tetranscripts:2.2.3--pyh7cba7a3_0' }"

    input:
    (meta_t, bam_t, bai_t): Tuple<Map, Path, Path?>
    (meta_c, bam_c, bai_c): Tuple<Map, Path, Path?>
    (meta_ggtf, g_gtf): Tuple<Map, Path>
    (meta_tegtf, te_gtf): Tuple<Map, Path>

    output:
    countTable= tuple(val(meta_t), files('*.cntTable'))
    log2fc= tuple(val(meta_t), files('*.R'))
    versions= file('versions.yml')
    DGE= tuple(val(meta_t), file('*.txt', optional: true))

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tetranscripts: \$(TEtranscripts --version)
END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta_c.id}"
    """
    echo $args
    
    touch ${prefix}.R
    touch ${prefix}.cntTable

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tetranscripts: \$(TEtranscripts --version)
    END_VERSIONS
    """
}
