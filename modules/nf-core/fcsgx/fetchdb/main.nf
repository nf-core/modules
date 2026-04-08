process FCSGX_FETCHDB {
    tag "$manifest.baseName"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.5--h9948957_0':
        'biocontainers/ncbi-fcs-gx:0.5.5--h9948957_0' }"

    input:
    val manifest // URL of manifest. Should not stage locally.

    output:
    path "$prefix"      , emit: database
    tuple val("${task.process}"), val('fcsgx'), eval("gx --help | sed '/build/!d; s/.*:v//; s/-.*//'"), emit: versions_fcsgx, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "gxdb_$manifest.baseName"
    """
    sync_files.py \\
        get \\
        --mft "${manifest.toUriString()}" \\
        --dir "$prefix"
    """

    stub:
    // def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "gxdb_$manifest.baseName"
    """
    touch ${prefix}
    """
}
