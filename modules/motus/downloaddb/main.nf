process MOTUS_DOWNLOADDB {
    label 'process_low'
    label 'error_ignore'

    conda (params.enable_conda ? "bioconda::motus=3.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:3.0.1--pyhdfd78af_0' :
        'quay.io/biocontainers/motus:3.0.1--pyhdfd78af_0' }"

    input:
    path motus_downloaddb

    output:
    path "db_mOTU"                                       , emit: db
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ## must copy file to working directory,
    ## otherwise the reference_db will be download to other than current directory
    cp $motus_downloaddb dwdDB.py
    python dwdDB.py \\
        -t $task.cpus
    rm dwdDB.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mOTUs: \$(grep motus db_mOTU/db_mOTU_versions | sed 's/motus\\t//g')
    END_VERSIONS
    """
}
