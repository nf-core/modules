process FGUMI_DEDUP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/99/99c92db2efcbcc4d20f2541060e09ed0f5a4338927f4a2ae3ab683de585efeb6/data'
        : 'community.wave.seqera.io/library/fgumi:0.7.0--d91f99b4cd4aae5a'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.metrics.txt"), emit: metrics
    tuple val(meta), path("*.family_size_histogram.txt"), emit: histogram
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_dedup"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    fgumi \\
        dedup \\
        --input ${bam} \\
        --output ${prefix}.bam \\
        --metrics ${prefix}.metrics.txt \\
        --family-size-histogram ${prefix}.family_size_histogram.txt \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_dedup"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.bam
    touch ${prefix}.metrics.txt
    touch ${prefix}.family_size_histogram.txt
    """
}
