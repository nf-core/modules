process KMCP_INDEX {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kmcp:0.9.4--h9ee0642_0'
        : 'quay.io/biocontainers/kmcp:0.9.4--h9ee0642_0'}"

    input:
    tuple val(meta), path(compute_dir)
    tuple val(meta2), path(taxdmp, stageAs: "taxdmp/*")
    tuple val(meta3), path(seq2taxidmap, stageAs: "taxdmp/*")

    output:
    tuple val(meta), path("${prefix}"), emit: kmcp
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('kmcp'), eval("kmcp version 2>&1 | sed 's/^.*kmcp v//'"), emit: versions_kmcp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    kmcp \\
        index \\
        --in-dir ${compute_dir} \\
        ${args} \\
        --threads ${task.cpus} \\
        --log ${prefix}.log \\
        --out-dir ${prefix}

    ## Optionally copy over taxonomy files if provided for downstream module
    ${taxdmp || seq2taxidmap ? "cp -r taxdmp/ ${prefix}/" : ''}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}.log

    ## Optionally copy over taxonomy files if provided for downstream module
    ${taxdmp || seq2taxidmap ? "cp -r taxdmp/ ${prefix}/" : ''}
    """
}
