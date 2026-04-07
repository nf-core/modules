process PIGZ_UNCOMPRESS {
    tag "$meta.id"
    label 'process_low'
    //stageInMode 'copy' // this directive can be set in case the original input should be kept

    conda "conda-forge::pigz"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.8':
        'biocontainers/pigz:2.8' }"

    input:
    tuple val(meta), path(zip)

    output:
    tuple val(meta), path("${uncompressed_filename}") , emit: file
    tuple val("${task.process}"), val('pigz'), eval('pigz --version 2>&1 | sed "s/^.*pigz[[:space:]]*//"'), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    uncompressed_filename = zip.toString() - '.gz'
    // calling pigz -f to make it follow symlinks
    """
    unpigz \\
        -p $task.cpus \\
        -fk \\
        $args \\
        ${zip}

    """

    stub:
    uncompressed_filename = zip.toString() - '.gz'
    """
    touch $uncompressed_filename
    """
}
