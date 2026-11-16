// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process RUNDBCAN_DATABASE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e8/e843a691d1c9be41394116f30ee3cf3f7bd042e0d59cb2f6a0a34b2b0a308e30/data'
        : 'community.wave.seqera.io/library/dbcan:5.2.9--fb040f67fa175a9b'}"

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
