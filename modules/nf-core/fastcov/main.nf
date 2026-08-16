process FASTCOV {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastcov:0.1.3--hdfd78af_0':
        'quay.io/biocontainers/fastcov:0.1.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(index)
    val(file_ext)

    output:
    tuple val(meta), path("${prefix}.${file_ext}")                                                                           , emit: coverage_plot
    tuple val(meta), path("log.txt")                                                                                         , emit: log
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('fastcov'), val('0.1.3')                                                               , emit: versions_fastcov, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    file_ext = file_ext.isEmpty() ? 'png' : file_ext
    """
    fastcov.py \\
        -o ${prefix}.${file_ext} \\
        ${args} \\
        ${bam} > log.txt 2>&1
    """

    stub:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    file_ext = file_ext.isEmpty() ? 'png' : file_ext

    """
    echo $args

    touch ${prefix}.${file_ext}
    touch log.txt
    """
}
