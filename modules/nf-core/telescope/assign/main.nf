process TELESCOPE_ASSIGN {
    tag "$meta_bam.id"
    label 'process_single'

    //conda "${moduleDir}/environment.yml" -- No conda at the moment
    container 'docker.io/hanalysis/telescope_1.0.3_clip'

    input:
    tuple val(meta_bam), path(bam)
    tuple val(meta_gtf), path(gtf)

    output:
    tuple val(meta_bam), path("*{updated,other}.bam"), emit: bam, optional: true // only for --updated_sam
    tuple val(meta_bam), path("*{updated,other}.sam"), emit: sam, optional: true // only for --updated_sam
    tuple val(meta_bam), path("*.tsv"), emit: tsv, optional: true // for when there's no alignments
    tuple val(meta_bam), path("*.log"), emit: log, optional: true
    tuple val("${task.process}"), val('telescope'), eval("telescope version | sed '1!d;s/.* //'"), emit: versions_telescope, topic: versions

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
        $args \\
        > telescope.log 2>&1

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

    """
}
