process RUNDBCAN_DATABASE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.2.6--pyhdfd78af_0' :
        'quay.io/biocontainers/dbcan:5.2.6--pyhdfd78af_0' }"

    output:
    path "dbcan_db", emit: dbcan_db
    tuple val("${task.process}"), val('rundbcan'), eval("run_dbcan version | sed 's/dbCAN version: //g'"), emit: versions_rundbcan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_dbcan database \\
        --db_dir dbcan_db \\
        --aws_s3
    """

    stub:
    """
    mkdir -p dbcan_db
    """
}
