process SINGLEM_DBDOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/singlem:0.20.3--pyhdfd78af_2'
        : 'biocontainers/singlem:0.20.3--pyhdfd78af_2'}"

    output:
    path("*.smpkg.zb")                                                                                                     , emit: singlem_database
    tuple val("${task.process}"), val('singlem'), eval('singlem --version'), topic: versions                                , emit: versions_singlem

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    singlem \\
        data \\
        --output-directory . \\
        ${args}
    """

    stub:
    """
    touch S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb
    """
}
