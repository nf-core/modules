process TELESCOPE_ASSIGN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ea/eaa68e02f66957e51ecba656dfeaa8576bada172780d95592105bc316a01a65c/data':
        'community.wave.seqera.io/library/telescope:1.0.3_fix--d176f12022b914cf' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*{updated,other}.bam"), emit: bam, optional: true // only for --updated_sam
    tuple val(meta), path("*{updated,other}.sam"), emit: sam, optional: true // only for --updated_sam
    tuple val(meta), path("*.tsv"), emit: tsv, optional: true // for when there's no alignments
    tuple val(meta), path("*.log"), emit: log, optional: true
    tuple val("${task.process}"), val('telescope'), eval("telescope --version | sed '1!d;s/.* //'"), emit: versions_telescope, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

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
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch ${prefix}-updated.bam
    touch ${prefix}-other.bam
    touch ${prefix}-updated.sam
    touch ${prefix}-other.sam
    touch ${prefix}-telescope_report.tsv
    """
}
