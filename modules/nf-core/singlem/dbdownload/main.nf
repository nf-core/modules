process SINGLEM_DBDOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'community.wave.seqera.io/library/singlem:0.20.3--52ee62e56540a33a' :
    'biocontainers/singlem:0.19.0--pyhdfd78af_0' }"

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
