process METHURATOR_GTESTIMATOR {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/methurator:2.1.1--pyhdfd78af_0'
        : 'biocontainers/methurator:2.1.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta3), path(fasta)

    output:
    tuple val(meta), path("methurator_summary.yml"), emit: summary_report
    tuple val("${task.process}"), val('methurator'), eval('methurator --version'), emit: versions_methurator, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if (params.rrbs) {
        args += '--rrbs'
    }
    """
    methurator gt-estimator \\
        ${bam} \\
        --fasta ${fasta} \\
        --minimum-coverage ${params.minimum_coverage} \\
        --t-max ${params.t_max} \\
        -@ ${task.cpus} \\
        --outdir . \\
        --compute_ci \\
        ${args}

    """

    stub:
    """
    touch methurator_summary.yml

    """
}
