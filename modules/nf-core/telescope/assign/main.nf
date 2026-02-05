process TELESCOPE_ASSIGN {
    tag "$meta_bam.id"
    label 'process_single'

    //conda "${moduleDir}/environment.yml" -- No conda at the moment
    container 'docker.io/hanalysis/telescope_1.0.3_clip'

    input:
    tuple val(meta_bam), path(bam)
    tuple val(meta_gtf), path(gtf)

    output:
    stdout
    tuple val(meta_bam), path("*{updated,other}.bam"), emit: bam, optional: true // only for --updated_sam
    tuple val(meta_bam), path("*{updated,other}.sam"), emit: sam, optional: true // only for --updated_sam
    tuple val(meta_bam), path("*.tsv"), emit: tsv, optional: true // for when there's no alignments
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta_bam.id}"

    """
    
    echo -n "" > telescope.log

    telescope \\
        assign \\
        $bam \\
        $gtf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telescope: \$(telescope --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta_bam.id}"

    """
    echo $args

    touch ${prefix}-updated.bam
    touch ${prefix}-other.bam
    touch ${prefix}-updated.sam
    touch ${prefix}-other.sam
    touch ${prefix}-telescope_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telescope: \$(telescope --version)
    END_VERSIONS
    """
}
