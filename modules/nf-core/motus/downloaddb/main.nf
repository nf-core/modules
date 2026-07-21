process MOTUS_DOWNLOADDB {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:3.1.0--pyhdfd78af_0':
        'quay.io/biocontainers/motus:3.1.0--pyhdfd78af_0' }"

    input:
    path motus_downloaddb_script

    output:
    path "db_mOTU/", emit: db
    // WARN: Version information not provided by tool on CLI.  Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('motus'), val("3.1.0"), topic: versions, emit: versions_motus

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def software = "${motus_downloaddb_script.simpleName}_copy.py"
    """
    ## must copy script file to working directory,
    ## otherwise the reference_db will be download to bin folder
    ## other than current directory
    cp ${motus_downloaddb_script} ${software}
    python ${software} \\
        ${args} \\
        -t ${task.cpus}
    """

    stub:
    """
    mkdir db_mOTU
    """

}
