process MSISENSOR2_MSI {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/msisensor2:0.1--hd03093a_0'
        : 'biocontainers/msisensor2:0.1--hd03093a_0'}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(intervals)
    path scan
    path models, stageAs: "models/*"

    output:
    tuple val(meta), path("${prefix}"),          emit: msi
    tuple val(meta), path("${prefix}_dis"),      emit: distribution
    tuple val(meta), path("${prefix}_germline"), emit: germline, optional: true
    tuple val(meta), path("${prefix}_somatic"),  emit: somatic
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def interval_command = intervals  ? "-e ${intervals}" : ""
    def model_cmd        = models     ? "-M models/" : ""
    def normal_bam_cmd   = normal_bam ? "-n ${normal_bam}" : ""
    def scan_cmd         = scan       ? "-d ${scan}" : ""
    def tumor_bam_cmd    = tumor_bam  ? "-t ${tumor_bam}" : ""
    """
    msisensor2 msi \\
        -b ${task.cpus} \\
        ${args} \\
        ${model_cmd} \\
        ${scan_cmd} \\
        ${interval_command} \\
        ${tumor_bam_cmd} \\
        ${normal_bam_cmd} \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor2: \$(echo \$(msisensor2 2> >(grep Version) | sed 's/Version: v//g'))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}
    touch ${prefix}_dis
    touch ${prefix}_somatic
    touch ${prefix}_germline

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor2: \$(echo \$(msisensor2 2> >(grep Version) | sed 's/Version: v//g'))
    END_VERSIONS
    """
}
