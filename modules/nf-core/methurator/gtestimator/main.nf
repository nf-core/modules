process METHURATOR_GTESTIMATOR {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/methurator:2.1.1--pyhdfd78af_0'
        : 'quay.io/biocontainers/methurator:2.1.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta)

    output:
    tuple val(meta), path("${prefix}.yml"), emit: summary_report
    tuple val("${task.process}"), val('methurator'), eval("methurator --version | sed 's/.* //'"), emit: versions_methurator, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    methurator gt-estimator \\
        ${bam} \\
        --fasta ${fasta} \\
        -@ ${task.cpus} \\
        --outdir . \\
        ${args}

    mv methurator_summary.yml ${prefix}.yml
    """

    stub:
    prefix = task.ext.prefix?: "${meta.id}"
    """
    touch ${prefix}.yml

    """
}
