process PMLST_PREPARE_PMLST_DB {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pmlst:2.0.3--hdfd78af_0':
        'biocontainers/pmlst:2.0.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(db_dir)

    output:
    tuple val(meta), path('pmlst_db'), emit: db
    tuple val("${task.process}"), val('pmlst'), eval('echo 2.0.3'), emit: versions_pmlst, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cp -a ${db_dir} prepared_db
    cd prepared_db
    python3 INSTALL.py kma_index
    cd ..
    mv prepared_db pmlst_db
    """

    stub:
    """
    mkdir -p pmlst_db
    touch pmlst_db/config
    """
}
