process MOTUS_DOWNLOADDB {
    label 'process_low'

    conda "bioconda::motus=3.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:3.0.3--pyhdfd78af_0':
        'quay.io/biocontainers/motus:3.0.3--pyhdfd78af_0' }"

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

    ## mOTUs version number is not available from command line.
    ## mOTUs save the version number in index database folder.
    ## mOTUs will check the database version is same version as exec version.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        motus: \$(grep motus db_mOTU/db_mOTU_versions | sed 's/motus\\t//g')
    END_VERSIONS
    """
}
