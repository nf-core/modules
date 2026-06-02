process IDR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/idr:2.0.4.2--py39hcbe4a3b_5' :
        'quay.io/biocontainers/idr:2.0.4.2--py39hcbe4a3b_5' }"

    input:
    tuple val(meta), path(peaks), val(peak_type)

    output:
    tuple val(meta), path("*idrValues.txt"), emit: idr
    tuple val(meta), path("*log.txt")      , emit: log
    tuple val(meta), path("*.png")         , emit: png
    tuple val("${task.process}"), val('idr'), eval("idr --version |& sed '1!d;s/^.*IDR //'"), emit: versions_idr, topic: versions
    tuple val("${task.process}"), val('python'), eval("python --version |& sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if (peaks.toList().size < 2) {
        log.error "[ERROR] idr needs at least two replicates only one provided."
    }
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    if (!peak_types.contains(peak_type)) {
        log.error "[ERROR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    idr \\
        --samples ${peaks} \\
        --input-file-type ${peak_type} \\
        --output-file ${prefix}.idrValues.txt \\
        --log-output-file ${prefix}.log.txt \\
        --plot \\
        ${args}
    """

    stub:
    if (peaks.toList().size < 2) {
        log.error "[ERROR] idr needs at least two replicates only one provided."
    }
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    if (!peak_types.contains(peak_type)) {
        log.error "[ERROR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.idrValues.txt
    touch ${prefix}.log.txt
    touch ${prefix}.png
    """
}
