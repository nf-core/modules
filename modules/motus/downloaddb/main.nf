process MOTUS_DOWNLOADDB {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::motus=3.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:3.0.1--pyhdfd78af_0':
        'quay.io/biocontainers/motus:3.0.1--pyhdfd78af_0' }"

    input:
    path motus_downloaddb_script

    output:
    path "db_mOTU/"                , emit: db
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def software = "${motus_downloaddb_script.simpleName}_copy.py"
    """
    ## must copy script file to working directory,
    ## otherwise the reference_db will be download to bin folder
    ## other than current directory
    cp $motus_downloaddb_script ${software}
    python ${software} \\
        $args \\
        -t $task.cpus
    ## clean up
    rm ${software}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mOTUs: \$(echo \$(motus -h 2>&1) | sed 's/^.*Version: //; s/References.*\$//')
    END_VERSIONS
    """
}
