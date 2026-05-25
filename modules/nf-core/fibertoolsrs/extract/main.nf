process FIBERTOOLSRS_EXTRACT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fibertools-rs:0.7.1--h3b373d1_0':
        'quay.io/biocontainers/fibertools-rs:0.7.1--h3b373d1_0' }"

    input:
    tuple val(meta), path(bam)
    val(extract_type)

    output:
    tuple val(meta), path("*.bed.gz"), emit: bed
    tuple val("${task.process}"), val('fibertools-rs'), eval("ft --version | sed 's/fibertools-rs v//;s/\\t.*//'"), topic: versions, emit: versions_fibertoolsrs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def allowed_types = ['m6a', 'cpg', 'msp', 'nuc', 'all']
    assert allowed_types.contains(extract_type) : "Invalid extract_type='${extract_type}'. Allowed: ${allowed_types.join(', ')}"

    """
    ft \\
        extract \\
        ${args} \\
        --threads ${task.cpus} \\
        ${bam} \\
        --${extract_type} \\
        ${prefix}.bed.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.bed.gz
    """
}
